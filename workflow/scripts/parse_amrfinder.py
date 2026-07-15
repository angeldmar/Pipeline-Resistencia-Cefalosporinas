"""Normalizacion de resultados de deteccion de AMR (AMRFinderPlus).

AMRFinderPlus reporta, por cada gen o mutacion de resistencia detectado en un
ensamblaje, una fila con muchas columnas propias de la herramienta. Este
script las reduce y renombra a un esquema propio, estable y facil de
combinar entre muestras (tabla LARGA: una fila por gen detectado, no una
columna por gen), tal como recomienda la seccion de integracion de
resultados del diseno del pipeline.

Subcomandos:

  parse   -> normaliza la salida de AMRFinderPlus de UNA muestra y escribe el
             resultado en results/tables/amr/{sample_id}.tsv (un archivo por
             muestra; si no se detecto ningun gen, se escribe igual un
             archivo con encabezado y cero filas, para dejar constancia de
             que la muestra si se proceso).

  combine -> junta todas las tablas por muestra en un unico listado largo:
             results/tables/amr_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re

import pandas as pd

AMR_LONG_TABLE_COLUMNS = [
    "sample_id",
    "gene_symbol",
    "gene_family",
    "allele",
    "sequence_name",
    "antimicrobial_class",
    "antimicrobial_subclass",
    "detection_method",
    "percent_identity",
    "percent_coverage",
    "meets_identity_coverage_threshold",
    "contig_id",
    "start",
    "stop",
]

# Heuristica para separar la "familia" de gen del sufijo de alelo (ej.
# "blaCTX-M-15" -> familia "blaCTX-M", alelo "blaCTX-M-15"). AMRFinderPlus no
# entrega familia y alelo como columnas separadas: su "Gene symbol" ya
# identifica el alelo especifico cuando el metodo de deteccion es ALLELEX o
# EXACTX. Esto es solo descriptivo (para la tabla de salida); la
# clasificacion BLEE/AmpC/carbapenemasa en classify_cephalosporin_genes.py
# compara directamente contra el simbolo completo, no contra esta familia.
ALLELE_SUFFIX_PATTERN = re.compile(r"[-_]\d+(\.\d+)?[A-Za-z]?$")


def derive_gene_family(gene_symbol: str) -> str:
    """Quita un sufijo de alelo numerico final (ej. "-15", "-1.1") del
    simbolo del gen, si lo tiene."""
    return ALLELE_SUFFIX_PATTERN.sub("", gene_symbol)


def load_amrfinder_table(amrfinder_output_path: Path) -> pd.DataFrame:
    """Carga la tabla cruda que produce AMRFinderPlus."""
    return pd.read_csv(amrfinder_output_path, sep="\t")


def normalize_amrfinder_table(
    sample_id: str,
    raw_amrfinder_table: pd.DataFrame,
    minimum_identity: float,
    minimum_gene_coverage: float,
) -> pd.DataFrame:
    """Convierte la tabla cruda de AMRFinderPlus de una muestra en filas del
    esquema largo del pipeline, marcando ademas si cada deteccion cumple los
    umbrales minimos de identidad/cobertura de config.yaml (sin descartar
    las que no los cumplen: quedan igual en la tabla, solo marcadas)."""
    if raw_amrfinder_table.empty:
        return pd.DataFrame(columns=AMR_LONG_TABLE_COLUMNS)

    normalized_rows = []
    for _, raw_row in raw_amrfinder_table.iterrows():
        gene_symbol = str(raw_row["Gene symbol"])
        percent_identity = float(raw_row["% Identity to reference sequence"])
        percent_coverage = float(raw_row["% Coverage of reference sequence"])

        normalized_rows.append({
            "sample_id": sample_id,
            "gene_symbol": gene_symbol,
            "gene_family": derive_gene_family(gene_symbol),
            "allele": gene_symbol,
            "sequence_name": raw_row["Sequence name"],
            "antimicrobial_class": raw_row["Class"],
            "antimicrobial_subclass": raw_row["Subclass"],
            "detection_method": raw_row["Method"],
            "percent_identity": percent_identity,
            "percent_coverage": percent_coverage,
            "meets_identity_coverage_threshold": (
                percent_identity >= minimum_identity and percent_coverage >= minimum_gene_coverage
            ),
            "contig_id": raw_row["Contig id"],
            "start": int(raw_row["Start"]),
            "stop": int(raw_row["Stop"]),
        })

    return pd.DataFrame(normalized_rows, columns=AMR_LONG_TABLE_COLUMNS)


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
            f"No se encontraron tablas de AMR por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values(["sample_id", "gene_symbol"]).reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    amrfinder_output_path: Path,
    output_dir: Path,
    minimum_identity: float,
    minimum_gene_coverage: float,
) -> None:
    raw_amrfinder_table = load_amrfinder_table(amrfinder_output_path)
    normalized_table = normalize_amrfinder_table(
        sample_id, raw_amrfinder_table, minimum_identity, minimum_gene_coverage
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
        description="Normalizar y combinar resultados de deteccion de AMR (AMRFinderPlus)."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Normalizar la salida de AMRFinderPlus de una muestra"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("amrfinder_output", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/amr"),
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
        "--input-dir", type=Path, default=Path("results/tables/amr"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/amr_summary.tsv"),
        help="Ruta del listado combinado",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.amrfinder_output, args.output_dir,
            args.minimum_identity, args.minimum_gene_coverage,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
