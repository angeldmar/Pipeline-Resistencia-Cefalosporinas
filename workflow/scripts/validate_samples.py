"""Validacion de la tabla de muestras (config/samples.tsv) del pipeline de AMR en E. coli.

Este script se ejecuta antes de cualquier descarga o procesamiento. Su objetivo es
detectar errores de metadatos lo antes posible, para no descubrirlos horas despues
en medio de un ensamblaje fallido. Si encuentra problemas, los junta todos y los
reporta en un solo mensaje (en vez de detenerse en el primero) para que se puedan
corregir de una sola vez.
"""

from pathlib import Path
import argparse
import re

import pandas as pd

# Columnas que deben existir siempre en samples.tsv. Si falta alguna, el pipeline
# no puede continuar porque etapas posteriores (descarga, reportes, comparacion
# contra el estandar de referencia) dependen de ellas.
REQUIRED_COLUMNS = {
    "sample_id",
    "run_accession",
    "biosample",
    "sequencing_platform",
    "phenotype_cefotaxime",
    "phenotype_ceftriaxone",
    "phenotype_ceftazidime",
    "expected_genes",
    "data_source",
}

# Columnas de fenotipo, que solo pueden tomar los valores S (susceptible),
# I (intermedio), R (resistente) o NA (no disponible/no evaluado).
PHENOTYPE_COLUMNS = [
    "phenotype_cefotaxime",
    "phenotype_ceftriaxone",
    "phenotype_ceftazidime",
]
VALID_PHENOTYPE_VALUES = {"S", "I", "R", "NA"}

# Columnas opcionales: si un usuario ya tiene lecturas locales (en vez de
# descargarlas desde SRA/ENA), puede indicar su ruta aqui y el validador
# comprobara que esos archivos existen en disco.
OPTIONAL_LOCAL_READ_COLUMNS = ["local_fastq_r1", "local_fastq_r2"]

# Formatos de accesion aceptados. Las corridas de secuenciacion de SRA/ENA/DDBJ
# comienzan con SRR, ERR o DRR seguido de digitos; los BioSample con SAMN, SAMEA
# o SAMD seguido de digitos.
RUN_ACCESSION_PATTERN = re.compile(r"^(SRR|ERR|DRR)\d+$")
BIOSAMPLE_PATTERN = re.compile(r"^SAM(N|EA|D)\d+$")


def load_samples_table(samples_path: Path) -> pd.DataFrame:
    """Lee samples.tsv como texto plano (dtype=str) para no perder ceros a la
    izquierda ni convertir accesiones en numeros, y normaliza celdas vacias a "NA"."""
    return pd.read_csv(samples_path, sep="\t", dtype=str).fillna("NA")


def check_required_columns(samples_table: pd.DataFrame) -> list[str]:
    """Verifica que no falte ninguna columna obligatoria."""
    missing_columns = REQUIRED_COLUMNS.difference(samples_table.columns)
    if missing_columns:
        return [f"Faltan columnas obligatorias: {sorted(missing_columns)}"]
    return []


def check_duplicate_sample_ids(samples_table: pd.DataFrame) -> list[str]:
    """Verifica que sample_id sea unico y estable (no se usa el nombre de
    archivo como identificador, por lo que no debe haber duplicados aqui)."""
    duplicated_mask = samples_table["sample_id"].duplicated()
    if duplicated_mask.any():
        duplicated_ids = samples_table.loc[duplicated_mask, "sample_id"].tolist()
        return [f"Identificadores de muestra duplicados: {duplicated_ids}"]
    return []


def check_valid_phenotypes(samples_table: pd.DataFrame) -> list[str]:
    """Verifica que las columnas de fenotipo solo contengan S, I, R o NA."""
    errors = []
    for phenotype_column in PHENOTYPE_COLUMNS:
        invalid_values = set(samples_table[phenotype_column]) - VALID_PHENOTYPE_VALUES
        if invalid_values:
            errors.append(
                f"Valores invalidos en {phenotype_column}: {sorted(invalid_values)}"
            )
    return errors


def check_accession_formats(samples_table: pd.DataFrame) -> list[str]:
    """Verifica que run_accession y biosample tengan un formato reconocible,
    para detectar errores de captura (typos, columnas desplazadas, etc.)
    antes de intentar descargar datos con una accesion invalida."""
    errors = []
    for _, row in samples_table.iterrows():
        sample_id = row["sample_id"]
        if not RUN_ACCESSION_PATTERN.match(row["run_accession"]):
            errors.append(
                f"{sample_id}: run_accession con formato invalido "
                f"({row['run_accession']!r}); se esperaba SRR/ERR/DRR + digitos"
            )
        if not BIOSAMPLE_PATTERN.match(row["biosample"]):
            errors.append(
                f"{sample_id}: biosample con formato invalido "
                f"({row['biosample']!r}); se esperaba SAMN/SAMEA/SAMD + digitos"
            )
    return errors


def check_documented_source(samples_table: pd.DataFrame) -> list[str]:
    """Verifica que cada muestra tenga una fuente de datos documentada
    (repositorio, estudio o articulo de origen), para mantener trazabilidad."""
    missing_source_mask = samples_table["data_source"].isin(["NA", ""])
    if missing_source_mask.any():
        samples_without_source = samples_table.loc[missing_source_mask, "sample_id"].tolist()
        return [f"Muestras sin fuente de datos documentada: {samples_without_source}"]
    return []


def check_local_read_files_exist(samples_table: pd.DataFrame) -> list[str]:
    """Si la tabla incluye columnas de lecturas locales (local_fastq_r1/r2),
    verifica que esos archivos realmente existan en disco. Estas columnas son
    opcionales: solo se usan cuando una muestra ya tiene FASTQ locales en vez
    de descargarse desde un repositorio publico."""
    available_columns = [
        column for column in OPTIONAL_LOCAL_READ_COLUMNS if column in samples_table.columns
    ]
    if not available_columns:
        return []

    errors = []
    for _, row in samples_table.iterrows():
        sample_id = row["sample_id"]
        for column in available_columns:
            file_path = row[column]
            if file_path in ("NA", ""):
                continue
            if not Path(file_path).is_file():
                errors.append(
                    f"{sample_id}: archivo indicado en {column} no existe ({file_path})"
                )
    return errors


def validate_samples(samples_path: Path) -> pd.DataFrame:
    """Ejecuta todas las verificaciones sobre samples.tsv y devuelve la tabla
    validada. Si se encuentra al menos un problema, lanza ValueError con la
    lista completa de errores encontrados."""
    samples_table = load_samples_table(samples_path)

    # La verificacion de columnas obligatorias debe pasar primero: si faltan
    # columnas, el resto de las verificaciones no tendria sentido (fallarian
    # con un KeyError en vez de un mensaje claro).
    column_errors = check_required_columns(samples_table)
    if column_errors:
        raise ValueError("\n".join(column_errors))

    all_errors: list[str] = []
    all_errors += check_duplicate_sample_ids(samples_table)
    all_errors += check_valid_phenotypes(samples_table)
    all_errors += check_accession_formats(samples_table)
    all_errors += check_documented_source(samples_table)
    all_errors += check_local_read_files_exist(samples_table)

    if all_errors:
        error_report = "\n".join(f"  - {error}" for error in all_errors)
        raise ValueError(f"La validacion de samples.tsv encontro errores:\n{error_report}")

    return samples_table


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validar la tabla de muestras del pipeline antes de descargar o procesar datos."
    )
    parser.add_argument("samples", type=Path, help="Ruta a config/samples.tsv")
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Ruta donde se escribira la tabla validada (TSV)",
    )
    args = parser.parse_args()

    validated_samples_table = validate_samples(args.samples)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    validated_samples_table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
