"""Descarga de lecturas crudas desde SRA para el pipeline de AMR en E. coli.

Este script controla (no reemplaza) las herramientas oficiales de descarga:
usa `fasterq-dump` (SRA Toolkit) via subprocess para obtener las lecturas
paired-end de cada muestra, las comprime a gzip, calcula su hash SHA-256 y
deja un manifiesto de descarga con la trazabilidad completa (muestra,
accesion, fuente, fecha, tamano, hash y estado de la descarga).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

# Las lecturas paired-end siempre tienen dos direcciones: forward (R1) y
# reverse (R2). El pipeline no acepta datos single-end.
READ_DIRECTIONS = ("R1", "R2")

# Etiqueta fija del repositorio de origen, para dejar constancia en el
# manifiesto de que las lecturas provienen de SRA via fasterq-dump.
REPOSITORY_LABEL = "NCBI SRA (fasterq-dump)"

# Columnas del manifiesto de descarga (una fila por archivo R1/R2).
MANIFEST_COLUMNS = [
    "sample_id",
    "run_accession",
    "read_direction",
    "repository",
    "data_source",
    "download_date",
    "file_path",
    "file_size_bytes",
    "sha256",
    "download_status",
    "error_message",
]


def run_command(command: list[str]) -> None:
    """Ejecuta un comando externo y lanza un error legible si falla.

    Python solo orquesta el proceso (via subprocess); la descarga real la
    hace la herramienta oficial (fasterq-dump), no un cliente HTTP casero.
    """
    try:
        completed_process = subprocess.run(command, check=False, capture_output=True, text=True)
    except FileNotFoundError as error:
        raise RuntimeError(
            f"No se encontro el ejecutable '{command[0]}'. "
            "Verifica que el ambiente conda con sra-tools este activo."
        ) from error

    if completed_process.returncode != 0:
        raise RuntimeError(
            f"Error ejecutando '{' '.join(command)}':\n{completed_process.stderr}"
        )


def download_sra_reads(run_accession: str, download_threads: int, temp_dir: Path) -> tuple[Path, Path]:
    """Descarga las lecturas paired-end de una accesion SRA con fasterq-dump.

    Devuelve las rutas a los dos archivos FASTQ sin comprimir (_1 y _2).
    Lanza RuntimeError si no se obtuvieron ambos archivos, ya que el pipeline
    requiere datos paired-end y no debe continuar silenciosamente con datos
    incompletos.
    """
    run_command([
        "fasterq-dump",
        run_accession,
        "--split-files",
        "--outdir", str(temp_dir),
        "--threads", str(download_threads),
    ])

    forward_read_path = temp_dir / f"{run_accession}_1.fastq"
    reverse_read_path = temp_dir / f"{run_accession}_2.fastq"

    if not forward_read_path.is_file() or not reverse_read_path.is_file():
        raise RuntimeError(
            f"No se obtuvieron ambos archivos paired-end para {run_accession} "
            f"(se esperaba {forward_read_path.name} y {reverse_read_path.name})"
        )
    return forward_read_path, reverse_read_path


def compress_fastq(source_path: Path, destination_path: Path) -> None:
    """Comprime un archivo FASTQ sin comprimir a formato gzip."""
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    with open(source_path, "rb") as source_file, gzip.open(destination_path, "wb") as destination_file:
        shutil.copyfileobj(source_file, destination_file)


def compute_sha256(file_path: Path) -> str:
    """Calcula el hash SHA-256 de un archivo leyendolo en bloques, para no
    cargar en memoria archivos FASTQ comprimidos que pueden pesar varios GB."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(8192), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def _build_manifest_row(
    sample_id: str,
    run_accession: str,
    read_direction: str,
    data_source: str,
    download_date: str,
    file_path: Path | None,
    download_status: str,
    error_message: str,
) -> dict:
    """Arma una fila del manifiesto con un formato consistente, tanto para
    descargas exitosas como fallidas."""
    return {
        "sample_id": sample_id,
        "run_accession": run_accession,
        "read_direction": read_direction,
        "repository": REPOSITORY_LABEL,
        "data_source": data_source,
        "download_date": download_date,
        "file_path": str(file_path) if file_path is not None else "",
        "file_size_bytes": file_path.stat().st_size if file_path is not None else 0,
        "sha256": compute_sha256(file_path) if file_path is not None else "",
        "download_status": download_status,
        "error_message": error_message,
    }


def download_sample(
    sample_id: str,
    run_accession: str,
    data_source: str,
    raw_data_dir: Path,
    download_threads: int,
) -> list[dict]:
    """Descarga, comprime y verifica las lecturas de una muestra.

    Devuelve dos filas de manifiesto (R1 y R2). Si la descarga falla en
    cualquier punto, ambas filas quedan marcadas como FAILED con el motivo
    del error, en vez de detener la descarga del resto de las muestras.
    """
    download_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    final_read_paths = {
        "R1": raw_data_dir / f"{sample_id}_R1.fastq.gz",
        "R2": raw_data_dir / f"{sample_id}_R2.fastq.gz",
    }

    try:
        with tempfile.TemporaryDirectory(prefix=f"{run_accession}_") as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            forward_read_path, reverse_read_path = download_sra_reads(
                run_accession, download_threads, temp_dir
            )
            compress_fastq(forward_read_path, final_read_paths["R1"])
            compress_fastq(reverse_read_path, final_read_paths["R2"])

        return [
            _build_manifest_row(
                sample_id, run_accession, read_direction, data_source, download_date,
                file_path=final_read_paths[read_direction],
                download_status="SUCCESS",
                error_message="",
            )
            for read_direction in READ_DIRECTIONS
        ]

    except (RuntimeError, OSError) as error:
        return [
            _build_manifest_row(
                sample_id, run_accession, read_direction, data_source, download_date,
                file_path=None,
                download_status="FAILED",
                error_message=str(error),
            )
            for read_direction in READ_DIRECTIONS
        ]


def update_manifest(manifest_path: Path, new_rows: list[dict]) -> pd.DataFrame:
    """Combina las filas nuevas con el manifiesto existente.

    Si una muestra ya tenia filas en el manifiesto (por ejemplo, de un
    intento de descarga anterior), esas filas se reemplazan por las nuevas
    en vez de acumularse como duplicados.
    """
    new_rows_table = pd.DataFrame(new_rows, columns=MANIFEST_COLUMNS)

    if manifest_path.is_file():
        existing_table = pd.read_csv(manifest_path, sep="\t", dtype=str)
        reprocessed_sample_ids = set(new_rows_table["sample_id"])
        existing_table = existing_table[~existing_table["sample_id"].isin(reprocessed_sample_ids)]
        combined_table = pd.concat([existing_table, new_rows_table], ignore_index=True)
    else:
        combined_table = new_rows_table

    return combined_table.sort_values(["sample_id", "read_direction"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Descargar lecturas paired-end desde SRA para las muestras del pipeline."
    )
    parser.add_argument("samples", type=Path, help="Ruta a samples.tsv (ya validado)")
    parser.add_argument(
        "--raw-dir", type=Path, default=Path("data/raw"),
        help="Carpeta de salida para los FASTQ comprimidos",
    )
    parser.add_argument(
        "--manifest", type=Path, default=Path("data/metadata/download_manifest.tsv"),
        help="Ruta del manifiesto de descarga",
    )
    parser.add_argument("--threads", type=int, default=4, help="Hilos usados por fasterq-dump")
    parser.add_argument(
        "--sample-id", type=str, default=None,
        help="Procesar unicamente esta muestra (por defecto: todas las de samples.tsv)",
    )
    args = parser.parse_args()

    samples_table = pd.read_csv(args.samples, sep="\t", dtype=str)
    if args.sample_id is not None:
        samples_table = samples_table[samples_table["sample_id"] == args.sample_id]
        if samples_table.empty:
            raise ValueError(f"No se encontro la muestra '{args.sample_id}' en {args.samples}")

    args.raw_dir.mkdir(parents=True, exist_ok=True)

    all_new_rows: list[dict] = []
    for _, sample_row in samples_table.iterrows():
        print(f"Descargando {sample_row['sample_id']} ({sample_row['run_accession']})...")
        all_new_rows.extend(
            download_sample(
                sample_id=sample_row["sample_id"],
                run_accession=sample_row["run_accession"],
                data_source=sample_row["data_source"],
                raw_data_dir=args.raw_dir,
                download_threads=args.threads,
            )
        )

    updated_manifest = update_manifest(args.manifest, all_new_rows)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    updated_manifest.to_csv(args.manifest, sep="\t", index=False)

    processed_sample_ids = set(samples_table["sample_id"])
    failed_mask = (
        updated_manifest["sample_id"].isin(processed_sample_ids)
        & (updated_manifest["download_status"] == "FAILED")
    )
    failed_sample_ids = sorted(updated_manifest.loc[failed_mask, "sample_id"].unique())
    if failed_sample_ids:
        print(f"Fallaron {len(failed_sample_ids)} muestra(s): {failed_sample_ids}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
