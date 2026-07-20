"""Extraccion de metricas de identificacion taxonomica desde reportes de Kraken2.

Kraken2 clasifica cada lectura contra una base de datos de referencia y
produce un reporte jerarquico (una fila por taxon, con el porcentaje de
lecturas asignadas a el y a su descendencia). Este script reduce ese reporte
a lo que el pipeline necesita para confirmar que una muestra es realmente
Escherichia coli: el taxon predominante, el % de lecturas asignado a E. coli,
y el % asignado a otras especies (posible contaminacion).

REGLA DE REVISION MANUAL (no automatizar la exclusion): Shigella es tan
proxima genomicamente a E. coli que ambos generos historicamente se
consideraron la misma especie; Kraken2 puede clasificar una fraccion de las
lecturas de una muestra de E. coli como Shigella (o viceversa) sin que eso
signifique contaminacion real. Por eso las lecturas asignadas a Shigella NO
se suman al porcentaje de "otras especies" que dispara un FAIL, y en cambio
la muestra se marca con requires_manual_review=True para que un analista
revise el caso manualmente, en vez de excluirla o aceptarla automaticamente.

La revision se dispara si el % de Shigella supera
shigella_review_threshold_percentage (config: taxonomy, por defecto 0.1%),
no con cualquier valor mayor a cero: en datos reales, Kraken2 casi siempre
asigna un rastro minimo de lecturas a generos cercanos por ambiguedad de
k-mers, sin que eso sea una senal real de mezcla de especies -- un umbral
en cero convertiria la revision manual en una alerta que se dispara casi
siempre, restandole utilidad como senal.

Subcomandos:

  parse   -> lee el report.tsv de Kraken2 de UNA muestra y escribe su fila de
             metricas en results/tables/taxonomy/{sample_id}.tsv.

  combine -> junta todas las tablas por muestra en
             results/tables/taxonomy_summary.tsv, y ademas escribe
             results/tables/taxonomy_manual_review.tsv con las muestras que
             tienen presencia de Shigella y requieren revision manual.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

# Columnas del reporte estandar de Kraken2 (sin encabezado en el archivo original).
KRAKEN2_REPORT_COLUMNS = [
    "percentage",
    "reads_in_clade",
    "reads_direct",
    "rank_code",
    "taxid",
    "name",
]

# Codigo de rango taxonomico que identifica una fila a nivel de especie.
SPECIES_RANK_CODE = "S"

TARGET_SPECIES_NAME = "Escherichia coli"
CLOSE_RELATIVE_GENUS_NAME = "Shigella"

TAXONOMY_SUMMARY_COLUMNS = [
    "sample_id",
    "predominant_taxon",
    "ecoli_percentage",
    "shigella_percentage",
    "other_contaminant_percentage",
    "requires_manual_review",
    "taxonomy_status",
]


def load_kraken2_report(kraken2_report_path: Path) -> pd.DataFrame:
    """Carga el reporte de Kraken2 y limpia los nombres de taxon (Kraken2 los
    indenta con espacios para reflejar la profundidad en el arbol taxonomico)."""
    report_table = pd.read_csv(
        kraken2_report_path, sep="\t", header=None, names=KRAKEN2_REPORT_COLUMNS
    )
    report_table["name"] = report_table["name"].str.strip()
    return report_table


def classify_taxonomy(
    ecoli_percentage: float,
    other_contaminant_percentage: float,
    minimum_ecoli_percentage: float,
    warning_ecoli_percentage: float,
    maximum_contaminant_percentage: float,
) -> str:
    """Clasifica la identificacion taxonomica en PASS, WARNING o FAIL.

    Nota: other_contaminant_percentage excluye deliberadamente a Shigella
    (ver docstring del modulo); esta funcion no sabe nada de Shigella, solo
    recibe el porcentaje de contaminantes ya depurado de ese caso especial.
    """
    if ecoli_percentage >= minimum_ecoli_percentage and other_contaminant_percentage < maximum_contaminant_percentage:
        return "PASS"
    if ecoli_percentage >= warning_ecoli_percentage:
        return "WARNING"
    return "FAIL"


def extract_taxonomy_metrics(
    sample_id: str,
    kraken2_report: pd.DataFrame,
    minimum_ecoli_percentage: float,
    warning_ecoli_percentage: float,
    maximum_contaminant_percentage: float,
    shigella_review_threshold_percentage: float = 0.1,
) -> dict:
    """Reduce el reporte completo de Kraken2 a las metricas de interes del
    pipeline, aplicando la regla especial de Shigella descrita arriba."""
    species_level_rows = kraken2_report.loc[kraken2_report["rank_code"] == SPECIES_RANK_CODE]

    ecoli_rows = species_level_rows.loc[species_level_rows["name"] == TARGET_SPECIES_NAME]
    ecoli_percentage = float(ecoli_rows["percentage"].iloc[0]) if not ecoli_rows.empty else 0.0

    shigella_rows = species_level_rows.loc[
        species_level_rows["name"].str.startswith(CLOSE_RELATIVE_GENUS_NAME)
    ]
    shigella_percentage = round(float(shigella_rows["percentage"].sum()), 2)
    # Un umbral minimo (no cero) evita que ruido tipico de clasificacion
    # (ambiguedad de k-mers entre generos cercanos, presente en casi toda
    # corrida real) dispare la revision manual de forma constante; ver
    # docstring del modulo.
    requires_manual_review = shigella_percentage > shigella_review_threshold_percentage

    other_species_rows = species_level_rows.loc[
        (species_level_rows["name"] != TARGET_SPECIES_NAME)
        & ~species_level_rows["name"].str.startswith(CLOSE_RELATIVE_GENUS_NAME)
    ]
    other_contaminant_percentage = round(float(other_species_rows["percentage"].sum()), 2)

    if species_level_rows.empty:
        predominant_taxon = "unclassified"
    else:
        predominant_taxon = species_level_rows.loc[
            species_level_rows["percentage"].idxmax(), "name"
        ]

    taxonomy_status = classify_taxonomy(
        ecoli_percentage, other_contaminant_percentage,
        minimum_ecoli_percentage, warning_ecoli_percentage, maximum_contaminant_percentage,
    )

    return {
        "sample_id": sample_id,
        "predominant_taxon": predominant_taxon,
        "ecoli_percentage": ecoli_percentage,
        "shigella_percentage": shigella_percentage,
        "other_contaminant_percentage": other_contaminant_percentage,
        "requires_manual_review": requires_manual_review,
        "taxonomy_status": taxonomy_status,
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de metricas de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=TAXONOMY_SUMMARY_COLUMNS).to_csv(
        output_path, sep="\t", index=False
    )
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de taxonomia por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def build_manual_review_registry(combined_table: pd.DataFrame) -> pd.DataFrame:
    """Extrae las muestras con presencia de Shigella para dejar constancia
    explicita de que requieren revision manual (no se auto-excluyen ni se
    auto-aprueban solo por esta senal)."""
    flagged_samples = combined_table.loc[combined_table["requires_manual_review"]].copy()
    return flagged_samples[
        ["sample_id", "predominant_taxon", "ecoli_percentage", "shigella_percentage", "taxonomy_status"]
    ].reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    kraken2_report_path: Path,
    output_dir: Path,
    minimum_ecoli_percentage: float,
    warning_ecoli_percentage: float,
    maximum_contaminant_percentage: float,
    shigella_review_threshold_percentage: float,
) -> None:
    kraken2_report = load_kraken2_report(kraken2_report_path)
    sample_metrics = extract_taxonomy_metrics(
        sample_id, kraken2_report,
        minimum_ecoli_percentage, warning_ecoli_percentage, maximum_contaminant_percentage,
        shigella_review_threshold_percentage,
    )
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"Metricas de {sample_id} escritas en {output_path} "
        f"(E. coli: {sample_metrics['ecoli_percentage']}%, "
        f"Shigella: {sample_metrics['shigella_percentage']}%, "
        f"estado: {sample_metrics['taxonomy_status']}, "
        f"revision manual: {sample_metrics['requires_manual_review']})"
    )


def run_combine_command(per_sample_dir: Path, summary_output_path: Path, manual_review_output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(summary_output_path, sep="\t", index=False)

    manual_review_registry = build_manual_review_registry(combined_table)
    manual_review_output_path.parent.mkdir(parents=True, exist_ok=True)
    manual_review_registry.to_csv(manual_review_output_path, sep="\t", index=False)

    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {summary_output_path}")
    print(
        f"Registro de revision manual ({len(manual_review_registry)} muestra(s) con Shigella) "
        f"escrito en {manual_review_output_path}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extraer y combinar metricas de identificacion taxonomica de reportes de Kraken2."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Extraer las metricas de una muestra desde su reporte de Kraken2"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("kraken2_report", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/taxonomy"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )
    parse_command.add_argument(
        "--minimum-ecoli-percentage", type=float, default=90.0,
        help="%% minimo asignado a E. coli para PASS (config: taxonomy.minimum_ecoli_percentage)",
    )
    parse_command.add_argument(
        "--warning-ecoli-percentage", type=float, default=70.0,
        help="%% minimo asignado a E. coli para WARNING en vez de FAIL (config: taxonomy.warning_ecoli_percentage)",
    )
    parse_command.add_argument(
        "--maximum-contaminant-percentage", type=float, default=5.0,
        help="%% maximo de otras especies (sin contar Shigella) para PASS (config: taxonomy.maximum_contaminant_percentage)",
    )
    parse_command.add_argument(
        "--shigella-review-threshold-percentage", type=float, default=0.1,
        help="%% minimo de Shigella para marcar revision manual, por encima de ruido tipico de "
             "clasificacion (config: taxonomy.shigella_review_threshold_percentage)",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un resumen y un registro de revision manual"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/taxonomy"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/taxonomy_summary.tsv"),
        help="Ruta de la tabla de resumen combinada",
    )
    combine_command.add_argument(
        "--manual-review-output", type=Path, default=Path("results/tables/taxonomy_manual_review.tsv"),
        help="Ruta del registro de muestras con presencia de Shigella que requieren revision manual",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.kraken2_report, args.output_dir,
            args.minimum_ecoli_percentage, args.warning_ecoli_percentage, args.maximum_contaminant_percentage,
            args.shigella_review_threshold_percentage,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output, args.manual_review_output)


if __name__ == "__main__":
    main()
