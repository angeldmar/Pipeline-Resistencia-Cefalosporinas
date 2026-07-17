"""Normalizacion de resultados del segundo motor de deteccion de AMR (ABricate).

ABricate es una herramienta de deteccion de genes de resistencia independiente
de AMRFinderPlus, usada en este pipeline exclusivamente para calcular
concordancia analitica entre dos motores distintos (ver
compare_amr_engines.py) -- no reemplaza a AMRFinderPlus como fuente principal
de deteccion, que sigue siendo la usada para comparar contra el estandar de
referencia fenotipico (compare_to_reference.py).

Este script reduce la salida cruda de ABricate al mismo esquema de columnas
que ya usa parse_amrfinder.py (tabla LARGA: una fila por gen detectado), para
que ambos motores queden en un formato directamente comparable.

Subcomandos:

  parse   -> normaliza la salida de ABricate de UNA muestra (puede incluir
             resultados de varias bases de datos, ej. CARD + ResFinder, ya
             concatenados) y escribe results/tables/abricate/{sample_id}.tsv

  combine -> junta todas las tablas por muestra en
             results/tables/abricate_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re

import pandas as pd

# Misma heuristica de familia de gen que parse_amrfinder.py (se duplica a
# proposito: cada script de este pipeline es autocontenido). Es especialmente
# importante aqui: CARD/ResFinder y el catalogo de referencia de NCBI
# (AMRFinderPlus) no siempre nombran el mismo gen de forma identica, asi que
# la comparacion entre motores se hace a nivel de familia, no de alelo exacto.
ALLELE_SUFFIX_PATTERN = re.compile(r"[-_]\d+(\.\d+)?[A-Za-z]?$")


def derive_gene_family(gene_symbol: str) -> str:
    return ALLELE_SUFFIX_PATTERN.sub("", gene_symbol)


ABRICATE_LONG_TABLE_COLUMNS = [
    "sample_id",
    "gene_symbol",
    "gene_family",
    "database",
    "product",
    "resistance_classes",
    "percent_identity",
    "percent_coverage",
    "meets_identity_coverage_threshold",
    "contig_id",
    "start",
    "stop",
]


def load_abricate_table(abricate_output_path: Path) -> pd.DataFrame:
    """Carga la tabla cruda que produce ABricate para UNA base de datos. Su
    encabezado real empieza con "#FILE"; se renombra a "FILE" para no
    arrastrar el simbolo "#"."""
    raw_table = pd.read_csv(abricate_output_path, sep="\t")
    return raw_table.rename(columns={"#FILE": "FILE"})


def load_abricate_tables(abricate_output_paths: list[Path]) -> pd.DataFrame:
    """Carga y concatena la salida cruda de ABricate de VARIAS bases de datos
    (ej. CARD + ResFinder) de la misma muestra en una sola tabla."""
    tables = [load_abricate_table(path) for path in abricate_output_paths]
    return pd.concat(tables, ignore_index=True)


def normalize_abricate_table(
    sample_id: str,
    raw_abricate_table: pd.DataFrame,
    minimum_identity: float,
    minimum_gene_coverage: float,
) -> pd.DataFrame:
    """Convierte la tabla cruda de ABricate (una o mas bases de datos
    concatenadas) en filas del esquema largo del pipeline, marcando (sin
    descartar) las detecciones que no alcanzan los umbrales configurados."""
    if raw_abricate_table.empty:
        return pd.DataFrame(columns=ABRICATE_LONG_TABLE_COLUMNS)

    normalized_rows = []
    for _, raw_row in raw_abricate_table.iterrows():
        gene_symbol = str(raw_row["GENE"])
        percent_identity = float(raw_row["%IDENTITY"])
        percent_coverage = float(raw_row["%COVERAGE"])

        normalized_rows.append({
            "sample_id": sample_id,
            "gene_symbol": gene_symbol,
            "gene_family": derive_gene_family(gene_symbol),
            "database": raw_row["DATABASE"],
            "product": raw_row.get("PRODUCT", ""),
            "resistance_classes": raw_row.get("RESISTANCE", ""),
            "percent_identity": percent_identity,
            "percent_coverage": percent_coverage,
            "meets_identity_coverage_threshold": (
                percent_identity >= minimum_identity and percent_coverage >= minimum_gene_coverage
            ),
            "contig_id": raw_row["SEQUENCE"],
            "start": int(raw_row["START"]),
            "stop": int(raw_row["END"]),
        })

    return pd.DataFrame(normalized_rows, columns=ABRICATE_LONG_TABLE_COLUMNS)


def write_per_sample_table(sample_id: str, normalized_table: pd.DataFrame, output_dir: Path) -> Path:
    """Escribe las filas de una muestra en su propio archivo TSV (incluso si
    quedaron 0 filas, para distinguir "sin genes detectados" de "no procesada")."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_id}.tsv"
    normalized_table.to_csv(output_path, sep="\t", index=False)
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en un unico listado largo, ordenado por muestra y gen."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de ABricate por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values(["sample_id", "gene_symbol"]).reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    abricate_output_paths: list[Path],
    output_dir: Path,
    minimum_identity: float,
    minimum_gene_coverage: float,
) -> None:
    raw_abricate_table = load_abricate_tables(abricate_output_paths)
    normalized_table = normalize_abricate_table(
        sample_id, raw_abricate_table, minimum_identity, minimum_gene_coverage
    )
    output_path = write_per_sample_table(sample_id, normalized_table, output_dir)
    print(f"{sample_id}: {len(normalized_table)} gen(es) de resistencia normalizados en {output_path}")


def run_combine_command(per_sample_dir: Path, output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(output_path, sep="\t", index=False)
    print(f"Listado combinado de {len(combined_table)} deteccion(es) escrito en {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Normalizar y combinar resultados del segundo motor de deteccion de AMR (ABricate)."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Normalizar la salida de ABricate de una muestra"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument(
        "abricate_outputs", type=Path, nargs="+",
        help="Uno o mas archivos crudos de ABricate (uno por base de datos) de la misma muestra",
    )
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/abricate"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )
    parse_command.add_argument(
        "--minimum-identity", type=float, default=90.0,
        help="%% identidad minima para marcar una deteccion como confiable (config: amr.minimum_identity)",
    )
    parse_command.add_argument(
        "--minimum-gene-coverage", type=float, default=80.0,
        help="%% cobertura minima para marcar una deteccion como confiable (config: amr.minimum_gene_coverage)",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un unico listado largo"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/abricate"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/abricate_summary.tsv"),
        help="Ruta del listado combinado",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.abricate_outputs, args.output_dir,
            args.minimum_identity, args.minimum_gene_coverage,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
