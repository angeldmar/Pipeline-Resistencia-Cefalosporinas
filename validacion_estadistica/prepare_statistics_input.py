"""Agregacion de resultados de la validacion estadistica en una sola tabla.

Distinto de workflow/scripts/merge_results.py (tabla maestra de la
herramienta operativa, con columnas de fastp/taxonomia/desempeno que esta
seccion no produce): este script es propio de validacion_estadistica/, junta
solo lo que run_validation_batch.py genero por muestra (comparacion con el
estandar de referencia, CheckM, QUAST, MLST) y le agrega la categoria de
control (ESBL/AmpC/NEGATIVO/CONTROL LIMITE) desde
muestras_validacion_completo.csv. Ningun calculo estadistico se hace aqui
-- eso es exclusivo del notebook de R (ver notebooks/validacion_estadistica.ipynb).

Uso:
    python prepare_statistics_input.py
Salida:
    resultados/tables/validation_summary.tsv (una fila por muestra)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

VALIDATION_DIR = Path(__file__).resolve().parent
RESULTS_DIR = VALIDATION_DIR / "resultados"
TABLES_DIR = RESULTS_DIR / "tables"
SAMPLES_CSV = VALIDATION_DIR / "muestras" / "muestras_validacion_completo.csv"
OUTPUT_PATH = TABLES_DIR / "validation_summary.tsv"


def concat_per_sample_tsvs(directory: Path) -> pd.DataFrame:
    tsv_files = sorted(directory.glob("*.tsv"))
    if not tsv_files:
        raise FileNotFoundError(f"No se encontraron tablas por muestra en {directory}")
    return pd.concat((pd.read_csv(tsv_file, sep="\t", dtype=str) for tsv_file in tsv_files), ignore_index=True)


def load_samples_metadata() -> pd.DataFrame:
    samples = pd.read_csv(SAMPLES_CSV, dtype=str)
    samples["sample_id"] = samples["Assembly Accession"].str.replace(".", "_", regex=False)
    return samples[
        ["sample_id", "Assembly Accession", "Clasificación", "Tipo de control", "Resultado conocido"]
    ].rename(columns={
        "Assembly Accession": "assembly_accession",
        "Clasificación": "category",
        "Tipo de control": "control_type",
        "Resultado conocido": "known_result",
    })


def main() -> None:
    reference_comparison = concat_per_sample_tsvs(TABLES_DIR / "reference_comparison")
    checkm = concat_per_sample_tsvs(TABLES_DIR / "checkm")
    quast = concat_per_sample_tsvs(TABLES_DIR / "quast")
    mlst = concat_per_sample_tsvs(TABLES_DIR / "mlst")
    samples_metadata = load_samples_metadata()

    summary = samples_metadata.merge(reference_comparison, on="sample_id", how="left")
    summary = summary.merge(
        checkm[["sample_id", "completeness_percent", "contamination_percent", "completeness_status"]],
        on="sample_id", how="left",
    )
    summary = summary.merge(
        quast[["sample_id", "contigs", "n50", "total_length", "assembly_status"]],
        on="sample_id", how="left",
    )
    summary = summary.merge(
        mlst[["sample_id", "sequence_type", "sequence_type_status"]],
        on="sample_id", how="left",
    )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT_PATH, sep="\t", index=False)

    print(f"{len(summary)} muestra(s) agregadas, escritas en {OUTPUT_PATH}")
    print(f"Distribucion confusion_category: {summary['confusion_category'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
