"""Combina los registros de desempeno (tiempo/CPU/RAM) de todas las
muestras y modulos en una sola tabla, y deriva resumenes utiles.

Cada ejecucion de run_with_timing.py deja un archivo
results/tables/performance/{sample_id}_{module}.tsv con una fila. Este
script los junta en:

  results/tables/performance_summary.tsv
      Listado largo, una fila por (muestra, modulo) -- tal como pide la
      seccion de medicion de desempeno del diseno del pipeline.

  results/tables/performance_by_sample.tsv
      Tiempo total y RAM maxima POR MUESTRA (sumando/maximizando sobre sus
      modulos). Estas son las dos columnas que se incorporan a la tabla
      maestra (merge_results.py).

  results/tables/performance_by_module.tsv
      Tiempo total y RAM maxima POR MODULO (sumando/maximizando sobre todas
      las muestras), para identificar el modulo mas costoso del pipeline.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd


def combine_per_execution_files(performance_dir: Path) -> pd.DataFrame:
    """Junta todos los archivos {sample_id}_{module}.tsv de un directorio en
    un unico listado largo, ordenado por muestra y modulo."""
    per_execution_files = sorted(performance_dir.glob("*.tsv"))
    if not per_execution_files:
        raise FileNotFoundError(
            f"No se encontraron registros de desempeno en {performance_dir}"
        )

    all_execution_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_execution_files]
    combined_table = pd.concat(all_execution_tables, ignore_index=True)
    return combined_table.sort_values(["sample_id", "module"]).reset_index(drop=True)


def summarize_by_sample(performance_table: pd.DataFrame) -> pd.DataFrame:
    """Tiempo TOTAL (suma de todos los modulos) y RAM MAXIMA (el pico entre
    todos los modulos, no la suma: los modulos no corren simultaneamente
    dentro de la misma muestra) por muestra."""
    return performance_table.groupby("sample_id").agg(
        total_elapsed_seconds=("elapsed_seconds", "sum"),
        peak_max_ram_gb=("max_ram_gb", "max"),
        modules_run=("module", "count"),
    ).reset_index()


def summarize_by_module(performance_table: pd.DataFrame) -> pd.DataFrame:
    """Tiempo total y RAM maxima por modulo, a traves de todas las muestras
    -- util para identificar el modulo mas costoso del pipeline y para el
    calculo posterior de coeficiente de variacion entre repeticiones."""
    return performance_table.groupby("module").agg(
        total_elapsed_seconds=("elapsed_seconds", "sum"),
        mean_elapsed_seconds=("elapsed_seconds", "mean"),
        peak_max_ram_gb=("max_ram_gb", "max"),
        samples_run=("sample_id", "nunique"),
    ).reset_index().sort_values("total_elapsed_seconds", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combinar los registros de desempeno del pipeline en resumenes por muestra y por modulo."
    )
    parser.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/performance"),
        help="Carpeta con los archivos {sample_id}_{module}.tsv individuales",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/tables/performance_summary.tsv"),
        help="Ruta del listado largo combinado (una fila por muestra+modulo)",
    )
    parser.add_argument(
        "--by-sample-output", type=Path, default=Path("results/tables/performance_by_sample.tsv"),
        help="Ruta del resumen de tiempo total / RAM maxima por muestra",
    )
    parser.add_argument(
        "--by-module-output", type=Path, default=Path("results/tables/performance_by_module.tsv"),
        help="Ruta del resumen de tiempo total / RAM maxima por modulo",
    )
    args = parser.parse_args()

    performance_table = combine_per_execution_files(args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    performance_table.to_csv(args.output, sep="\t", index=False)

    by_sample_summary = summarize_by_sample(performance_table)
    by_sample_summary.to_csv(args.by_sample_output, sep="\t", index=False)

    by_module_summary = summarize_by_module(performance_table)
    by_module_summary.to_csv(args.by_module_output, sep="\t", index=False)

    print(f"{len(performance_table)} registro(s) de desempeno combinados en {args.output}")
    print(f"Resumen por muestra escrito en {args.by_sample_output}")
    print(f"Resumen por modulo escrito en {args.by_module_output}")
    if not by_module_summary.empty:
        costliest_module = by_module_summary.iloc[0]
        print(
            f"Modulo mas costoso: {costliest_module['module']} "
            f"({costliest_module['total_elapsed_seconds']}s totales)"
        )


if __name__ == "__main__":
    main()
