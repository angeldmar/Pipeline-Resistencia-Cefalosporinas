"""Verificacion taxonomica (Kraken2) de las muestras de validacion estadistica.

Distinto de la regla kraken2 de la seccion operativa (workflow/rules/
taxonomy.smk), que clasifica LECTURAS crudas ya recortadas: esta seccion
solo tiene los genomas ya ensamblados (ver validacion_estadistica/README.md),
asi que Kraken2 se corre directamente sobre el ensamblaje (un fasta de
contigs en vez de un par de fastq). El formato del reporte de Kraken2 es el
mismo en ambos casos, asi que parse_kraken2.py (parse/combine) se reutiliza
sin cambios.

Motivacion: al revisar los 13 falsos negativos de reference_comparison,
varios mostraban %GC muy alejado del rango normal de E. coli (~50-51%) y
MLST sin ningun alelo resuelto -- fuerte indicio de genoma mal etiquetado
(no es realmente E. coli), no de una falla de deteccion de AMRFinderPlus.
Este script confirma o descarta esa sospecha con una herramienta de
identificacion taxonomica real, en vez de quedarse con el proxy GC+MLST.

Reanudable: igual criterio que run_validation_batch.py, si la tabla de
taxonomia de una muestra ya existe se salta.

Uso:
    python run_taxonomy_check.py [--workers N] [--threads-per-sample N]
Salida:
    resultados/tables/taxonomy/{sample_id}.tsv (una fila por muestra)
    resultados/tables/taxonomy_summary.tsv (combinado, via parse_kraken2.py combine)
    resultados/tables/taxonomy_manual_review.tsv (casos con Shigella)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
VALIDATION_DIR = Path(__file__).resolve().parent
GENOMES_DIR = VALIDATION_DIR / "muestras" / "genomas" / "ncbi_dataset" / "data"
RESULTS_DIR = VALIDATION_DIR / "resultados"
SAMPLES_CSV = VALIDATION_DIR / "muestras" / "muestras_validacion_completo.csv"

sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

SNAKEMAKE_CONDA_DIR = REPO_ROOT / ".snakemake" / "conda"


class MissingCondaEnvironmentError(RuntimeError):
    pass


def resolve_conda_env_bin(env_yaml_relative_path: str) -> Path:
    """Duplicado a proposito de run_validation_batch.py: cada script de esta
    seccion es autocontenido (ver su docstring)."""
    source_content = (REPO_ROOT / env_yaml_relative_path).read_text()
    if SNAKEMAKE_CONDA_DIR.is_dir():
        for marker in SNAKEMAKE_CONDA_DIR.glob("*.yaml"):
            if marker.read_text() == source_content:
                env_hash = marker.stem
                bin_dir = SNAKEMAKE_CONDA_DIR / env_hash / "bin"
                if (SNAKEMAKE_CONDA_DIR / f"{env_hash}.env_setup_done").is_file() and bin_dir.is_dir():
                    return bin_dir
    raise MissingCondaEnvironmentError(f"No se encontro un ambiente Conda ya creado para {env_yaml_relative_path}.")


def env_with_conda_bin(bin_dir: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    return env


def find_genome_fasta(accession: str) -> Path:
    candidates = list((GENOMES_DIR / accession).glob("*_genomic.fna"))
    if not candidates:
        raise FileNotFoundError(f"No se encontro el archivo genomico para {accession} en {GENOMES_DIR / accession}")
    return candidates[0]


def run_and_log(command: list[str], log_path: Path, env: dict[str, str] | None = None) -> None:
    with open(log_path, "a") as log_handle:
        log_handle.write(f"$ {' '.join(command)}\n")
        completed = subprocess.run(command, cwd=REPO_ROOT, stdout=log_handle, stderr=subprocess.STDOUT, env=env)
    if completed.returncode != 0:
        raise RuntimeError(f"'{' '.join(command)}' fallo con codigo {completed.returncode} (ver {log_path})")


def process_sample(sample_id: str, accession: str, threads: int, config: dict) -> tuple[str, str]:
    tables_dir = RESULTS_DIR / "tables" / "taxonomy"
    final_table = tables_dir / f"{sample_id}.tsv"
    if final_table.is_file():
        return sample_id, "ya procesada, se salta"

    log_path = RESULTS_DIR / "logs" / f"taxonomy_{sample_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("")

    try:
        genome_path = find_genome_fasta(accession)
        kraken2_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/kraken2.yaml"))

        kraken2_dir = RESULTS_DIR / "taxonomy" / sample_id
        kraken2_dir.mkdir(parents=True, exist_ok=True)
        report_path = kraken2_dir / "report.tsv"
        classification_path = kraken2_dir / "classification.tsv"

        # Sin --paired: se clasifica el ensamblaje (un fasta de contigs), no
        # un par de fastq de lecturas -- unica diferencia real frente a la
        # regla kraken2 de la seccion operativa (workflow/rules/taxonomy.smk).
        run_and_log([
            "kraken2", "--db", config["paths"]["kraken_database"],
            str(genome_path), "--threads", str(threads),
            "--report", str(report_path), "--output", str(classification_path),
        ], log_path, env=kraken2_env)

        run_and_log([
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_kraken2.py"), "parse", sample_id, str(report_path),
            "--output-dir", str(tables_dir),
            "--minimum-ecoli-percentage", str(config["taxonomy"]["minimum_ecoli_percentage"]),
            "--warning-ecoli-percentage", str(config["taxonomy"]["warning_ecoli_percentage"]),
            "--maximum-contaminant-percentage", str(config["taxonomy"]["maximum_contaminant_percentage"]),
            "--shigella-review-threshold-percentage", str(config["taxonomy"]["shigella_review_threshold_percentage"]),
        ], log_path)

        return sample_id, "OK"
    except Exception as error:  # noqa: BLE001 -- se registra en el log de la muestra, no se deja morir el hilo en silencio
        with open(log_path, "a") as log_handle:
            log_handle.write(f"ERROR: {error}\n")
        return sample_id, f"FALLO: {error}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Verificar taxonomicamente (Kraken2) las muestras de validacion estadistica.")
    parser.add_argument("--workers", type=int, default=2, help="Muestras procesadas en paralelo")
    parser.add_argument("--threads-per-sample", type=int, default=4, help="Hilos de Kraken2 por muestra")
    parser.add_argument("--limit", type=int, default=None, help="Procesar solo las primeras N muestras (para pruebas)")
    args = parser.parse_args()

    with open(REPO_ROOT / "config" / "config.yaml") as config_file:
        config = yaml.safe_load(config_file)

    samples_df = pd.read_csv(SAMPLES_CSV)
    if args.limit:
        samples_df = samples_df.head(args.limit)

    jobs = []
    for _, row in samples_df.iterrows():
        accession = row["Assembly Accession"]
        sample_id = accession.replace(".", "_")
        jobs.append((sample_id, accession))

    print(f"Verificando taxonomia de {len(jobs)} muestras con {args.workers} en paralelo...", flush=True)
    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_sample, sample_id, accession, args.threads_per_sample, config): sample_id
            for sample_id, accession in jobs
        }
        for future in as_completed(futures):
            sample_id, status = future.result()
            completed_count += 1
            print(f"[{completed_count}/{len(jobs)}] {sample_id}: {status}", flush=True)

    tables_dir = RESULTS_DIR / "tables" / "taxonomy"
    subprocess.run([
        sys.executable, str(REPO_ROOT / "workflow/scripts/parse_kraken2.py"), "combine",
        "--input-dir", str(tables_dir),
        "--output", str(RESULTS_DIR / "tables" / "taxonomy_summary.tsv"),
        "--manual-review-output", str(RESULTS_DIR / "tables" / "taxonomy_manual_review.tsv"),
    ], cwd=REPO_ROOT, check=True)


if __name__ == "__main__":
    main()
