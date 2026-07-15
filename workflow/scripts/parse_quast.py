"""Extraccion de metricas de ensamblaje desde los reportes de QUAST.

QUAST evalua un ensamblaje y entrega un report.tsv con muchas metricas (una
por fila). Este script se queda solo con las que el pipeline necesita para
decidir si un ensamblaje es utilizable (numero de contigs, contig mas largo,
longitud total, contenido GC y N50), y clasifica cada ensamblaje en
PASS/WARNING/FAIL segun los umbrales de config.yaml.

Igual que parse_fastp.py, tiene dos subcomandos:

  parse   -> lee el report.tsv de UNA muestra y escribe su fila de metricas
             en results/tables/quast/{sample_id}.tsv (un archivo por muestra,
             para evitar que corridas paralelas de Snakemake escriban a la
             vez sobre un mismo archivo compartido).

  combine -> junta todas las tablas por muestra en un unico resumen:
             results/tables/quast_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

# Nombres exactos de las filas de interes en report.tsv de QUAST. QUAST
# tambien incluye filas calificadas como "# contigs (>= 1000 bp)" o
# "Total length (>= 5000 bp)"; estas claves sin calificador identifican las
# filas que resumen el ensamblaje completo, sin filtrar por longitud.
QUAST_ROW_CONTIGS = "# contigs"
QUAST_ROW_LARGEST_CONTIG = "Largest contig"
QUAST_ROW_TOTAL_LENGTH = "Total length"
QUAST_ROW_GC_CONTENT = "GC (%)"
QUAST_ROW_N50 = "N50"

# Columnas de la tabla de resumen de QUAST, en el orden en que se reportan.
QUAST_SUMMARY_COLUMNS = [
    "sample_id",
    "contigs",
    "largest_contig",
    "total_length",
    "gc_content_percent",
    "n50",
    "assembly_status",
]


def load_quast_report(quast_report_path: Path) -> pd.Series:
    """Carga report.tsv de QUAST como una Serie (nombre de metrica -> valor).

    QUAST genera una columna por ensamblaje evaluado; como este pipeline
    corre QUAST una muestra a la vez, solo existe una columna de datos, que
    es la que se devuelve.
    """
    report_table = pd.read_csv(quast_report_path, sep="\t", index_col=0)
    single_assembly_column = report_table.columns[0]
    return report_table[single_assembly_column]


def classify_assembly(
    contigs: int,
    total_length: int,
    n50: int,
    maximum_contigs: int,
    minimum_total_length: int,
    maximum_total_length: int,
    n50_warning_threshold: int,
) -> str:
    """Clasifica un ensamblaje en PASS, WARNING o FAIL.

    - Demasiados contigs, o una longitud total fuera del rango esperado para
      un genoma de E. coli, indican un ensamblaje poco confiable -> FAIL.
    - Un N50 bajo indica un ensamblaje mas fragmentado de lo ideal, pero
      todavia utilizable -> WARNING.
    - En cualquier otro caso, el ensamblaje se considera aceptable -> PASS.
    """
    if contigs > maximum_contigs:
        return "FAIL"
    if not (minimum_total_length <= total_length <= maximum_total_length):
        return "FAIL"
    if n50 < n50_warning_threshold:
        return "WARNING"
    return "PASS"


def extract_quast_metrics(
    sample_id: str,
    quast_report: pd.Series,
    maximum_contigs: int,
    minimum_total_length: int,
    maximum_total_length: int,
    n50_warning_threshold: int,
) -> dict:
    """Reduce el reporte completo de QUAST a las metricas que le interesan al
    pipeline, y les asigna un estado PASS/WARNING/FAIL."""
    contig_count = int(quast_report[QUAST_ROW_CONTIGS])
    largest_contig_length = int(quast_report[QUAST_ROW_LARGEST_CONTIG])
    total_length = int(quast_report[QUAST_ROW_TOTAL_LENGTH])
    gc_content_percent = float(quast_report[QUAST_ROW_GC_CONTENT])
    n50 = int(quast_report[QUAST_ROW_N50])

    assembly_status = classify_assembly(
        contig_count, total_length, n50,
        maximum_contigs, minimum_total_length, maximum_total_length, n50_warning_threshold,
    )

    return {
        "sample_id": sample_id,
        "contigs": contig_count,
        "largest_contig": largest_contig_length,
        "total_length": total_length,
        "gc_content_percent": gc_content_percent,
        "n50": n50,
        "assembly_status": assembly_status,
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de metricas de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=QUAST_SUMMARY_COLUMNS).to_csv(
        output_path, sep="\t", index=False
    )
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de QUAST por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    quast_report_path: Path,
    output_dir: Path,
    maximum_contigs: int,
    minimum_total_length: int,
    maximum_total_length: int,
    n50_warning_threshold: int,
) -> None:
    quast_report = load_quast_report(quast_report_path)
    sample_metrics = extract_quast_metrics(
        sample_id, quast_report,
        maximum_contigs, minimum_total_length, maximum_total_length, n50_warning_threshold,
    )
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"Metricas de {sample_id} escritas en {output_path} "
        f"(contigs: {sample_metrics['contigs']}, N50: {sample_metrics['n50']}, "
        f"estado: {sample_metrics['assembly_status']})"
    )


def run_combine_command(per_sample_dir: Path, output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(output_path, sep="\t", index=False)
    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extraer y combinar metricas de ensamblaje de reportes de QUAST."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Extraer las metricas de una muestra desde su report.tsv de QUAST"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("quast_report", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/quast"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )
    parse_command.add_argument(
        "--maximum-contigs", type=int, default=500,
        help="Mas contigs que esto -> FAIL (config: assembly.maximum_contigs)",
    )
    parse_command.add_argument(
        "--minimum-total-length", type=int, default=4_000_000,
        help="Longitud total minima esperada -> config: assembly.minimum_total_length",
    )
    parse_command.add_argument(
        "--maximum-total-length", type=int, default=6_500_000,
        help="Longitud total maxima esperada -> config: assembly.maximum_total_length",
    )
    parse_command.add_argument(
        "--n50-warning-threshold", type=int, default=20_000,
        help="N50 por debajo de esto -> WARNING (config: assembly.n50_warning_threshold)",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un unico resumen"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/quast"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/quast_summary.tsv"),
        help="Ruta de la tabla de resumen combinada",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.quast_report, args.output_dir,
            args.maximum_contigs, args.minimum_total_length,
            args.maximum_total_length, args.n50_warning_threshold,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
