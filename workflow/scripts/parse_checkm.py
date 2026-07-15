"""Extraccion de metricas de completitud y contaminacion desde CheckM.

CheckM estima, a partir de genes marcadores de copia unica especificos del
linaje taxonomico, que tan completo esta un ensamblaje y cuanta secuencia
"de mas" (posible contaminacion o mezcla de cepas) contiene. Este script lee
la tabla que produce CheckM (`checkm lineage_wf --tab_table`) y clasifica
cada muestra en PASS/FAIL segun los umbrales de config.yaml.

Una muestra que falla NO se descarta en silencio: sigue apareciendo en el
resumen combinado con estado FAIL, y ademas queda registrada en un archivo
aparte de exclusiones (checkm_exclusions.tsv) para que quede explicito por
que se excluiria del analisis principal.

Subcomandos:

  parse   -> lee el report.tsv de CheckM de UNA muestra y escribe su fila de
             metricas en results/tables/checkm/{sample_id}.tsv (un archivo
             por muestra, para evitar condiciones de carrera en Snakemake).

  combine -> junta todas las tablas por muestra en results/tables/checkm_summary.tsv
             y escribe results/tables/checkm_exclusions.tsv con las muestras FAIL.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

# Nombres de columnas tal como los reporta "checkm lineage_wf --tab_table".
CHECKM_COLUMN_BIN_ID = "Bin Id"
CHECKM_COLUMN_COMPLETENESS = "Completeness"
CHECKM_COLUMN_CONTAMINATION = "Contamination"

CHECKM_SUMMARY_COLUMNS = [
    "sample_id",
    "completeness_percent",
    "contamination_percent",
    "completeness_status",
]

EXCLUSION_REASON = "completitud/contaminacion fuera de umbral (CheckM)"


def load_checkm_report(checkm_report_path: Path) -> pd.DataFrame:
    """Carga la tabla tabulada de CheckM (una fila por bin/muestra)."""
    return pd.read_csv(checkm_report_path, sep="\t")


def classify_completeness(
    completeness_percent: float,
    contamination_percent: float,
    minimum_completeness: float,
    maximum_contamination: float,
) -> str:
    """PASS solo si la completitud alcanza el minimo Y la contaminacion no
    supera el maximo; cualquier otro caso es FAIL."""
    if completeness_percent >= minimum_completeness and contamination_percent < maximum_contamination:
        return "PASS"
    return "FAIL"


def extract_checkm_metrics(
    sample_id: str,
    checkm_report: pd.DataFrame,
    minimum_completeness: float,
    maximum_contamination: float,
) -> dict:
    """Busca la fila de la muestra en la tabla de CheckM (por Bin Id) y
    calcula su estado de completitud/contaminacion."""
    sample_row = checkm_report.loc[checkm_report[CHECKM_COLUMN_BIN_ID] == sample_id]
    if sample_row.empty:
        raise ValueError(
            f"No se encontro la muestra '{sample_id}' en el reporte de CheckM "
            f"(Bin Id disponibles: {checkm_report[CHECKM_COLUMN_BIN_ID].tolist()})"
        )

    completeness_percent = float(sample_row.iloc[0][CHECKM_COLUMN_COMPLETENESS])
    contamination_percent = float(sample_row.iloc[0][CHECKM_COLUMN_CONTAMINATION])
    completeness_status = classify_completeness(
        completeness_percent, contamination_percent, minimum_completeness, maximum_contamination
    )

    return {
        "sample_id": sample_id,
        "completeness_percent": completeness_percent,
        "contamination_percent": contamination_percent,
        "completeness_status": completeness_status,
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de metricas de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=CHECKM_SUMMARY_COLUMNS).to_csv(
        output_path, sep="\t", index=False
    )
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de CheckM por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def build_exclusion_registry(combined_table: pd.DataFrame) -> pd.DataFrame:
    """Extrae las muestras FAIL del resumen combinado para dejar un registro
    explicito de exclusiones, en vez de simplemente descartarlas."""
    failed_samples = combined_table.loc[combined_table["completeness_status"] == "FAIL"].copy()
    failed_samples["reason"] = EXCLUSION_REASON
    return failed_samples[
        ["sample_id", "completeness_percent", "contamination_percent", "reason"]
    ].reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    checkm_report_path: Path,
    output_dir: Path,
    minimum_completeness: float,
    maximum_contamination: float,
) -> None:
    checkm_report = load_checkm_report(checkm_report_path)
    sample_metrics = extract_checkm_metrics(
        sample_id, checkm_report, minimum_completeness, maximum_contamination
    )
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"Metricas de {sample_id} escritas en {output_path} "
        f"(completitud: {sample_metrics['completeness_percent']}%, "
        f"contaminacion: {sample_metrics['contamination_percent']}%, "
        f"estado: {sample_metrics['completeness_status']})"
    )


def run_combine_command(per_sample_dir: Path, summary_output_path: Path, exclusions_output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(summary_output_path, sep="\t", index=False)

    exclusion_registry = build_exclusion_registry(combined_table)
    exclusions_output_path.parent.mkdir(parents=True, exist_ok=True)
    exclusion_registry.to_csv(exclusions_output_path, sep="\t", index=False)

    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {summary_output_path}")
    print(f"Registro de exclusiones ({len(exclusion_registry)} muestra(s) FAIL) escrito en {exclusions_output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extraer y combinar metricas de completitud/contaminacion de CheckM."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Extraer las metricas de una muestra desde el reporte tabular de CheckM"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("checkm_report", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/checkm"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )
    parse_command.add_argument(
        "--minimum-completeness", type=float, default=95.0,
        help="Completitud minima (%%) para PASS (config: assembly.minimum_completeness)",
    )
    parse_command.add_argument(
        "--maximum-contamination", type=float, default=5.0,
        help="Contaminacion maxima (%%) para PASS (config: assembly.maximum_contamination)",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un resumen y un registro de exclusiones"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/checkm"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/checkm_summary.tsv"),
        help="Ruta de la tabla de resumen combinada",
    )
    combine_command.add_argument(
        "--exclusions-output", type=Path, default=Path("results/tables/checkm_exclusions.tsv"),
        help="Ruta del registro de muestras excluidas por completitud/contaminacion",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.checkm_report, args.output_dir,
            args.minimum_completeness, args.maximum_contamination,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output, args.exclusions_output)


if __name__ == "__main__":
    main()
