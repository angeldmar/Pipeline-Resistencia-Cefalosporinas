"""Extraccion de metricas de control de calidad desde los reportes JSON de fastp.

fastp entrega, por cada muestra, un reporte JSON con estadisticas de las
lecturas antes y despues del recorte/filtrado. Este script reduce ese JSON a
las metricas que el pipeline necesita para juzgar la calidad de una muestra,
estima la cobertura del genoma a partir de esas mismas lecturas (sin correr
ninguna herramienta adicional), y deja todo en tablas TSV faciles de leer y
comparar entre muestras.

La cobertura NO se usa para descartar muestras en silencio: toda muestra
recibe un estado PASS/WARNING/FAIL segun los umbrales de config.yaml
(quality.minimum_coverage / quality.warning_coverage), y ese estado queda
registrado en la misma tabla de salida para que la decision de excluir o no
una muestra sea explicita y revisable.

Este script tiene dos modos de uso (subcomandos):

  parse   -> lee el JSON de UNA muestra y escribe su fila de metricas en
             results/tables/fastp/{sample_id}.tsv (un archivo por muestra,
             para que corridas paralelas de Snakemake no escriban a la vez
             sobre el mismo archivo).

  combine -> junta todas las tablas por muestra en un unico resumen:
             results/tables/fastp_summary.tsv
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json

import pandas as pd

# Columnas de la tabla de resumen de fastp, en el orden en que se reportan.
FASTP_SUMMARY_COLUMNS = [
    "sample_id",
    "initial_reads",
    "retained_reads",
    "filtered_reads_percent",
    "q20_bases",
    "q30_bases",
    "gc_content_percent",
    "duplication_rate_percent",
    "mean_read_length",
    "estimated_coverage",
    "coverage_status",
]


def load_fastp_report(fastp_json_path: Path) -> dict:
    """Carga el reporte JSON generado por fastp."""
    with open(fastp_json_path) as json_file:
        return json.load(json_file)


def estimate_coverage(read_count: int, mean_read_length: float, genome_size: int = 5_000_000) -> float:
    """Estima la cobertura del genoma a partir del numero de lecturas retenidas
    y su longitud media: Cobertura = (lecturas x longitud media) / tamano del genoma.

    Es una aproximacion rapida (no reemplaza medir la cobertura real sobre el
    ensamblaje), util para detectar temprano una muestra con datos insuficientes.
    """
    if genome_size <= 0:
        raise ValueError("El tamano del genoma debe ser mayor que cero.")
    return (read_count * mean_read_length) / genome_size


def classify_coverage(estimated_coverage_value: float, minimum_coverage: float, warning_coverage: float) -> str:
    """Clasifica la cobertura estimada en PASS, WARNING o FAIL.

    Una muestra con cobertura baja no se descarta silenciosamente: queda
    marcada para que un analista decida si la excluye o la mantiene con
    advertencia (ver seccion de control de calidad del diseno del pipeline).
    """
    if estimated_coverage_value >= minimum_coverage:
        return "PASS"
    if estimated_coverage_value >= warning_coverage:
        return "WARNING"
    return "FAIL"


def extract_fastp_metrics(
    sample_id: str,
    fastp_report: dict,
    genome_size: int,
    minimum_coverage: float,
    warning_coverage: float,
) -> dict:
    """Reduce el reporte completo de fastp a las metricas que le interesan al
    pipeline: lecturas iniciales/retenidas, bases de alta calidad (Q20/Q30),
    contenido GC, duplicacion y cobertura estimada, todas medidas despues del
    filtrado (sobre las lecturas que realmente se usaran en el ensamblaje)."""
    reads_before_filtering = fastp_report["summary"]["before_filtering"]
    reads_after_filtering = fastp_report["summary"]["after_filtering"]

    initial_read_count = reads_before_filtering["total_reads"]
    retained_read_count = reads_after_filtering["total_reads"]

    if initial_read_count == 0:
        filtered_reads_percent = 0.0
    else:
        reads_removed = initial_read_count - retained_read_count
        filtered_reads_percent = round(100 * reads_removed / initial_read_count, 2)

    # fastp reporta "rate" y "gc_content" como fracciones (0-1); se convierten
    # a porcentaje para que sean directamente comparables entre muestras.
    duplication_rate_percent = round(fastp_report["duplication"]["rate"] * 100, 2)
    gc_content_percent = round(reads_after_filtering["gc_content"] * 100, 2)

    # Longitud media de lectura tras el filtrado: se deriva de bases/lecturas
    # totales (promedio de R1 y R2 combinados) en vez de asumir un valor fijo.
    if retained_read_count == 0:
        mean_read_length = 0.0
        estimated_coverage_value = 0.0
    else:
        mean_read_length = round(reads_after_filtering["total_bases"] / retained_read_count, 2)
        estimated_coverage_value = round(
            estimate_coverage(retained_read_count, mean_read_length, genome_size), 2
        )

    coverage_status = classify_coverage(estimated_coverage_value, minimum_coverage, warning_coverage)

    return {
        "sample_id": sample_id,
        "initial_reads": initial_read_count,
        "retained_reads": retained_read_count,
        "filtered_reads_percent": filtered_reads_percent,
        "q20_bases": reads_after_filtering["q20_bases"],
        "q30_bases": reads_after_filtering["q30_bases"],
        "gc_content_percent": gc_content_percent,
        "duplication_rate_percent": duplication_rate_percent,
        "mean_read_length": mean_read_length,
        "estimated_coverage": estimated_coverage_value,
        "coverage_status": coverage_status,
    }


def write_per_sample_table(sample_metrics: dict, output_dir: Path) -> Path:
    """Escribe la fila de metricas de una muestra en su propio archivo TSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{sample_metrics['sample_id']}.tsv"
    pd.DataFrame([sample_metrics], columns=FASTP_SUMMARY_COLUMNS).to_csv(
        output_path, sep="\t", index=False
    )
    return output_path


def combine_per_sample_tables(per_sample_dir: Path) -> pd.DataFrame:
    """Junta todas las tablas individuales (una por muestra) de un directorio
    en una sola tabla de resumen, ordenada por sample_id."""
    per_sample_files = sorted(per_sample_dir.glob("*.tsv"))
    if not per_sample_files:
        raise FileNotFoundError(
            f"No se encontraron tablas de fastp por muestra en {per_sample_dir}"
        )

    all_sample_tables = [pd.read_csv(file_path, sep="\t") for file_path in per_sample_files]
    combined_table = pd.concat(all_sample_tables, ignore_index=True)
    return combined_table.sort_values("sample_id").reset_index(drop=True)


def run_parse_command(
    sample_id: str,
    fastp_json_path: Path,
    output_dir: Path,
    genome_size: int,
    minimum_coverage: float,
    warning_coverage: float,
) -> None:
    fastp_report = load_fastp_report(fastp_json_path)
    sample_metrics = extract_fastp_metrics(
        sample_id, fastp_report, genome_size, minimum_coverage, warning_coverage
    )
    output_path = write_per_sample_table(sample_metrics, output_dir)
    print(
        f"Metricas de {sample_id} escritas en {output_path} "
        f"(cobertura estimada: {sample_metrics['estimated_coverage']}x, "
        f"estado: {sample_metrics['coverage_status']})"
    )


def run_combine_command(per_sample_dir: Path, output_path: Path) -> None:
    combined_table = combine_per_sample_tables(per_sample_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined_table.to_csv(output_path, sep="\t", index=False)
    print(f"Resumen combinado de {len(combined_table)} muestra(s) escrito en {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extraer y combinar metricas de calidad de reportes JSON de fastp."
    )
    subcommands = parser.add_subparsers(dest="command", required=True)

    parse_command = subcommands.add_parser(
        "parse", help="Extraer las metricas de una muestra desde su reporte JSON de fastp"
    )
    parse_command.add_argument("sample_id", type=str)
    parse_command.add_argument("fastp_json", type=Path)
    parse_command.add_argument(
        "--output-dir", type=Path, default=Path("results/tables/fastp"),
        help="Carpeta donde se escribe la tabla individual de la muestra",
    )
    parse_command.add_argument(
        "--genome-size", type=int, default=5_000_000,
        help="Tamano de genoma esperado (pb) para estimar la cobertura (config: quality.estimated_genome_size)",
    )
    parse_command.add_argument(
        "--minimum-coverage", type=float, default=30.0,
        help="Cobertura minima (X) para marcar la muestra como PASS (config: quality.minimum_coverage)",
    )
    parse_command.add_argument(
        "--warning-coverage", type=float, default=15.0,
        help="Cobertura minima (X) para marcar la muestra como WARNING en vez de FAIL (config: quality.warning_coverage)",
    )

    combine_command = subcommands.add_parser(
        "combine", help="Combinar las tablas individuales en un unico resumen"
    )
    combine_command.add_argument(
        "--input-dir", type=Path, default=Path("results/tables/fastp"),
        help="Carpeta con las tablas individuales por muestra",
    )
    combine_command.add_argument(
        "--output", type=Path, default=Path("results/tables/fastp_summary.tsv"),
        help="Ruta de la tabla de resumen combinada",
    )

    args = parser.parse_args()

    if args.command == "parse":
        run_parse_command(
            args.sample_id, args.fastp_json, args.output_dir,
            args.genome_size, args.minimum_coverage, args.warning_coverage,
        )
    elif args.command == "combine":
        run_combine_command(args.input_dir, args.output)


if __name__ == "__main__":
    main()
