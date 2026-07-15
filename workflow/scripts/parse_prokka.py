"""Organizacion de resultados de anotacion genomica (Prokka).

Prokka anota un ensamblaje y produce varios archivos (.gff, .gbk, .faa,
.ffn, .tsv, .txt, entre otros). Este script no vuelve a anotar nada: solo
lee dos de esos archivos para dejar un resumen numerico y trazable por
muestra:

  - {sample}.txt  -> resumen que el propio Prokka calcula (CDS, rRNA, tRNA, ...).
  - {sample}.tsv  -> tabla de features (una fila por gen anotado), usada aqui
                      unicamente para contar cuantos CDS quedaron anotados
                      como "hypothetical protein" (sin funcion conocida).

La version de Prokka se registra consultando "prokka --version" en el
momento de correr el script, para trazabilidad (no se asume una version fija).

Subcomandos:

  parse   -> procesa UNA muestra y escribe su fila de metricas en
             results/tables/annotation/{sample_id}.tsv.

  combine -> junta todas las tablas por muestra en
             results/tables/annotation_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import subprocess

import pandas as pd

ANNOTATION_SUMMARY_COLUMNS = [
    "sample_id",
    "cds_count",
    "rrna_count",
    "trna_count",
    "hypothetical_gene_count",
    "prokka_version",
]

# Nombre de la columna de producto y el valor exacto que Prokka usa para
# marcar un gen sin funcion conocida, en el archivo {sample}.tsv.
PROKKA_PRODUCT_COLUMN = "product"
PROKKA_FEATURE_TYPE_COLUMN = "ftype"
HYPOTHETICAL_PRODUCT_LABEL = "hypothetical protein"
CDS_FEATURE_TYPE = "CDS"


def parse_prokka_summary(summary_txt_path: Path) -> dict:
    """Lee el resumen que Prokka ya calcula ({sample}.txt), con lineas en
    formato "clave: valor" (ej. "CDS: 4823"), y devuelve los conteos de
    interes como enteros."""
    counts_by_key = {}
    with open(summary_txt_path) as summary_file:
        for line in summary_file:
            if ":" not in line:
                continue
            key, value = line.strip().split(":", maxsplit=1)
            counts_by_key[key.strip()] = value.strip()

    return {
        "cds_count": int(counts_by_key.get("CDS", 0)),
        "rrna_count": int(counts_by_key.get("rRNA", 0)),
        "trna_count": int(counts_by_key.get("tRNA", 0)),
    }


def count_hypothetical_genes(annotation_tsv_path: Path) -> int:
    """Cuenta cuantos CDS de la tabla de anotacion de Prokka quedaron
    marcados como "hypothetical protein" (sin funcion conocida asignada)."""
    annotation_table = pd.read_csv(annotation_tsv_path, sep="\t")
    is_cds = annotation_table[PROKKA_FEATURE_TYPE_COLUMN] == CDS_FEATURE_TYPE
    is_hypothetical = annotation_table[PROKKA_PRODUCT_COLUMN] == HYPOTHETICAL_PRODUCT_LABEL
    return int((is_cds & is_hypothetical).sum())


def get_prokka_version() -> str:
    """Consulta la version de Prokka instalada, para dejar trazabilidad de
    con que version se genero la anotacion. Si Prokka no esta disponible en
    el entorno actual, devuelve "unknown" en vez de fallar (la version es
    metadata util pero no bloqueante para el resto del pipeline)."""
    try:
        completed_process = subprocess.run(
            ["prokka", "--version"], capture_output=True, text=True, check=False
        )
    except FileNotFoundError:
        return "unknown"

    # Prokka imprime la version a stderr, con formato "prokka 1.14.6".
    version_output = (completed_process.stdout + completed_process.stderr).strip()
    return version_output if version_output else "unknown"


def extract_annotation_metrics(sample_id: str, summary_txt_path: Path, annotation_tsv_path: Path) -> dict:
    """Combina el resumen de Prokka, el conteo de genes hipoteticos y la
    version de la herramienta en una sola fila de metricas por muestra."""
    summary_counts = parse_prokka_summary(summary_txt_path)
    hypothetical_gene_count = count_hypothetical_genes(annotation_tsv_path)
    prokka_version = get_prokka_version()

    return {
        "sample_id": sample_id,
        "cds_count": summary_counts["cds_count"],
        "rrna_count": summary_counts["rrna_count"],
        "trna_count": summary_counts["trna_count"],
        "hypothetical_gene_count": hypothetical_gene_count,
        "prokka_version": prokka_version,
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de metricas de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=ANNOTATION_SUMMARY_COLUMNS).to_csv(
        output_path, sep="\t", index=False
    )
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de anotacion por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def run_parse_command(sample_id: str, summary_txt_path: Path, annotation_tsv_path: Path, output_dir: Path) -> None:
    sample_metrics = extract_annotation_metrics(sample_id, summary_txt_path, annotation_tsv_path)
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"Metricas de {sample_id} escritas en {output_path} "
        f"(CDS: {sample_metrics['cds_count']}, rRNA: {sample_metrics['rrna_count']}, "
        f"tRNA: {sample_metrics['trna_count']}, hipoteticos: {sample_metrics['hypothetical_gene_count']})"
    )


def run_combine_command(per_sample_dir: Path, output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(output_path, sep="\t", index=False)
    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Organizar y resumir resultados de anotacion genomica (Prokka)."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Extraer las metricas de una muestra desde su resumen y tabla de Prokka"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("--summary-txt", type=Path, required=True, help="Archivo {sample}.txt de Prokka")
    parse_command.add_argument("--annotation-tsv", type=Path, required=True, help="Archivo {sample}.tsv de Prokka")
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/annotation"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un unico resumen"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/annotation"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/annotation_summary.tsv"),
        help="Ruta de la tabla de resumen combinada",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(args.sample_id, args.summary_txt, args.annotation_tsv, args.output_dir)
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
