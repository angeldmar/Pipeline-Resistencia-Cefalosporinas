"""Prepara el archivo de entrada para el analisis estadistico en R.

Junta la comparacion contra el estandar de referencia (TP/TN/FP/FN, ver
compare_to_reference.py) con las metricas de desempeno por corrida (ver
combine_performance.py) en un unico CSV limpio:

    results/statistics/validation_input.csv

con las columnas exactas que pide la seccion de estadistica del diseno del
pipeline: sample_id, reference_result, pipeline_result, run, elapsed_seconds,
max_ram_gb.

De aqui en adelante, TODO el calculo estadistico (sensibilidad, especificidad,
kappa, intervalos de confianza, coeficiente de variacion) es responsabilidad
exclusiva de R (ver workflow/scripts/run_statistics.R); este script solo
prepara y limpia los datos, no calcula ninguna estadistica.

Las muestras sin estandar de referencia documentado (reference_status =
"indeterminate", ver compare_to_reference.py) se excluyen: no hay nada
contra que compararlas, y R espera un factor binario positive/negative.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re

import pandas as pd

VALIDATION_INPUT_COLUMNS = [
    "sample_id",
    "reference_result",
    "pipeline_result",
    "run",
    "elapsed_seconds",
    "max_ram_gb",
]

# Misma convencion de nombre de corrida de reproducibilidad que
# assess_reproducibility.py ("EC001_run2" -> muestra base "EC001", corrida 2).
# Se duplica aqui a proposito: cada script de este pipeline es autocontenido.
REPLICATE_RUN_ID_PATTERN = re.compile(r"^(?P<base_sample_id>.+)_run(?P<run_number>\d+)$")


def parse_sample_and_run(sample_id: str) -> tuple[str, int]:
    """Separa un sample_id en (muestra_base, numero_de_corrida). Las
    muestras que no siguen la convencion "{base}_runN" (la gran mayoria,
    evaluadas una sola vez) se tratan como corrida 1 de si mismas."""
    match = REPLICATE_RUN_ID_PATTERN.match(sample_id)
    if match:
        return match.group("base_sample_id"), int(match.group("run_number"))
    return sample_id, 1


def build_validation_input(
    reference_comparison_table: pd.DataFrame,
    performance_table: pd.DataFrame,
) -> pd.DataFrame:
    """Combina la comparacion de referencia con el desempeno por corrida en
    el esquema que espera run_statistics.R."""
    scoreable_rows = reference_comparison_table.loc[
        reference_comparison_table["reference_status"] != "indeterminate"
    ].copy()

    parsed_ids = scoreable_rows["sample_id"].apply(parse_sample_and_run)
    scoreable_rows["base_sample_id"] = parsed_ids.apply(lambda pair: pair[0])
    scoreable_rows["run"] = parsed_ids.apply(lambda pair: pair[1])

    merged_table = scoreable_rows.merge(
        performance_table[["sample_id", "total_elapsed_seconds", "peak_max_ram_gb"]],
        on="sample_id", how="left",
    )

    validation_input = pd.DataFrame({
        "sample_id": merged_table["base_sample_id"],
        "reference_result": merged_table["reference_status"],
        "pipeline_result": merged_table["pipeline_status"],
        "run": merged_table["run"],
        "elapsed_seconds": merged_table["total_elapsed_seconds"],
        "max_ram_gb": merged_table["peak_max_ram_gb"],
    })

    return validation_input.sort_values(["sample_id", "run"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preparar validation_input.csv para el analisis estadistico en R."
    )
    parser.add_argument(
        "--reference-comparison", type=Path, default=Path("results/tables/reference_comparison.tsv"),
    )
    parser.add_argument(
        "--performance-by-sample", type=Path, default=Path("results/tables/performance_by_sample.tsv"),
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/statistics/validation_input.csv"),
    )
    args = parser.parse_args()

    reference_comparison_table = pd.read_csv(args.reference_comparison, sep="\t")
    performance_table = pd.read_csv(args.performance_by_sample, sep="\t")

    validation_input = build_validation_input(reference_comparison_table, performance_table)

    excluded_count = len(reference_comparison_table) - len(validation_input)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    validation_input.to_csv(args.output, index=False)

    print(f"{len(validation_input)} fila(s) escritas en {args.output}")
    if excluded_count:
        print(
            f"{excluded_count} fila(s) excluidas por no tener estandar de referencia documentado "
            "(reference_status='indeterminate')"
        )


if __name__ == "__main__":
    main()
