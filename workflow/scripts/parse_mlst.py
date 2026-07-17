"""Normalizacion de resultados de tipificacion de secuencia multilocus (MLST).

MLST da contexto epidemiologico (linaje/secuencia tipo, ej. ST131) a los
genes de resistencia ya detectados: dos aislamientos con el mismo gen de
resistencia pero distinto ST pueden representar eventos de transmision
independientes, mientras que el mismo ST en varias muestras sugiere
diseminacion clonal de un linaje. Este modulo no interactua con la deteccion
de AMR; solo enriquece la tabla maestra con una columna adicional de contexto.

La herramienta `mlst` (T. Seemann) imprime, por defecto, una fila por archivo
de entrada SIN encabezado: FILE, SCHEME, ST, y luego una columna por gen del
esquema con su alelo (ej. "adk(6)"). El numero de columnas de alelo varia
segun el esquema, asi que este script no asume una cantidad fija: junta todas
las columnas de alelo restantes en un solo texto ("allele_profile").

Un ST puede quedar como "-" cuando la combinacion de alelos no coincide con
ningun tipo conocido en la base de datos (sequence_type_status="unresolved"),
o marcarse con "~"/"?" en un alelo individual cuando la coincidencia es
aproximada, no exacta (sequence_type_status="novel_allele"). Ninguno de estos
casos se descarta ni se trata como error: quedan visibles en la tabla, tal
como el resto de los modulos de este pipeline nunca ocultan una muestra
problematica.

Subcomandos:

  parse   -> procesa la salida de mlst de UNA muestra y escribe su fila en
             results/tables/mlst/{sample_id}.tsv

  combine -> junta todas las tablas por muestra en
             results/tables/mlst_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

MLST_SUMMARY_COLUMNS = ["sample_id", "scheme", "sequence_type", "sequence_type_status", "allele_profile"]

# Marcadores que "mlst" usa para senalar una coincidencia de alelo no exacta
# (aproximada/novedosa), en vez de un alelo identico a uno ya catalogado.
INEXACT_ALLELE_MARKERS = ("~", "?")

UNRESOLVED_SEQUENCE_TYPE_VALUE = "-"


def load_mlst_output(mlst_output_path: Path) -> pd.DataFrame:
    """Carga la salida cruda de mlst (una fila por muestra, sin encabezado)."""
    return pd.read_csv(mlst_output_path, sep="\t", header=None)


def classify_sequence_type_status(sequence_type: str, allele_calls: list[str]) -> str:
    """Determina que tan confiable es el tipo de secuencia asignado:
    "unresolved" si no hubo ST catalogado, "novel_allele" si algun alelo
    individual solo tuvo una coincidencia aproximada, "exact" en el resto."""
    if sequence_type == UNRESOLVED_SEQUENCE_TYPE_VALUE:
        return "unresolved"
    has_inexact_allele = any(
        marker in allele_call for allele_call in allele_calls for marker in INEXACT_ALLELE_MARKERS
    )
    return "novel_allele" if has_inexact_allele else "exact"


def normalize_mlst_output(sample_id: str, raw_mlst_table: pd.DataFrame) -> dict:
    """Reduce la fila cruda de mlst de una muestra al esquema del pipeline."""
    if raw_mlst_table.empty:
        # No deberia ocurrir en condiciones normales (mlst siempre produce
        # una fila por archivo de entrada, incluso con ST no resuelto), pero
        # se maneja sin fallar en vez de asumir que nunca pasara.
        return {
            "sample_id": sample_id,
            "scheme": "unknown",
            "sequence_type": UNRESOLVED_SEQUENCE_TYPE_VALUE,
            "sequence_type_status": "unresolved",
            "allele_profile": "none",
        }

    raw_row = raw_mlst_table.iloc[0]
    scheme = str(raw_row[1])
    sequence_type = str(raw_row[2])
    allele_calls = [str(value) for value in raw_row[3:].dropna().tolist()]

    return {
        "sample_id": sample_id,
        "scheme": scheme,
        "sequence_type": sequence_type,
        "sequence_type_status": classify_sequence_type_status(sequence_type, allele_calls),
        "allele_profile": "; ".join(allele_calls) if allele_calls else "none",
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=MLST_SUMMARY_COLUMNS).to_csv(output_path, sep="\t", index=False)
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(f"No se encontraron tablas de MLST por muestra en {per_sample_dir}")

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def run_parse_command(sample_id: str, mlst_output_path: Path, output_dir: Path) -> None:
    raw_mlst_table = load_mlst_output(mlst_output_path)
    sample_metrics = normalize_mlst_output(sample_id, raw_mlst_table)
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"{sample_id}: esquema {sample_metrics['scheme']}, ST {sample_metrics['sequence_type']} "
        f"({sample_metrics['sequence_type_status']}), escrito en {output_path}"
    )


def run_combine_command(per_sample_dir: Path, output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(output_path, sep="\t", index=False)
    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalizar y combinar resultados de tipificacion MLST."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser("parse", help="Normalizar la salida de mlst de una muestra")
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("mlst_output", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/mlst"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )

    combine_command = subcommands.add_parser("combine", help="Combinar las tablas individuales en un resumen")
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/mlst"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/mlst_summary.tsv"),
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(args.sample_id, args.mlst_output, args.output_dir)
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
