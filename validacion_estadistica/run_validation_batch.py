"""Corrida por lotes del pipeline (modo solo-ensamblaje) sobre las muestras
de validacion estadistica.

Distinto de webapp/pipeline_runner.py (una carga ad-hoc a la vez) y del
Snakefile principal (lote curado de config/samples.tsv, con lecturas
crudas): este script procesa MUCHAS muestras ya ensambladas en paralelo,
cada una con las mismas herramientas y rutas de archivo que usarian las
reglas equivalentes de Snakemake (QUAST, CheckM, AMRFinderPlus, ABricate,
MLST), escribiendo a validacion_estadistica/resultados/ en vez de results/
(nunca se mezcla esta validacion con el lote curado de la herramienta
operativa).

Reanudable: si la tabla individual final de una muestra ya existe, esa
muestra se salta -- pausar (Ctrl+C o matar el proceso) y volver a correr
este script mas tarde retoma sin repetir trabajo ya hecho.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

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
    """Igual que en webapp/pipeline_runner.py: ubica el bin/ del ambiente
    Conda ya creado por Snakemake para esa herramienta, comparando el
    contenido del yaml de origen contra los marcadores en .snakemake/conda/
    (duplicado a proposito, mismo criterio de scripts autocontenidos)."""
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


# Valores de "Resultado conocido" (Excel/CSV) que significan "sin gen
# esperado" (referencia negativa) o "sin dato de referencia utilizable"
# (controles limite, no se evaluan por genotipo). El resto de valores se
# usan como expected_genes tal cual, o se reduce a la familia (sin sufijo
# de alelo) cuando el valor es una variante ambigua ("CTX-M-type ESBL*").
NEGATIVE_RESULT_MARKERS = ("ausencia esperada",)
NO_REFERENCE_RESULT_MARKERS = ("no aplica",)


def derive_expected_genes(resultado_conocido: str) -> str:
    """Traduce la columna 'Resultado conocido' del CSV de muestras al mismo
    formato de expected_genes que usa compare_to_reference.py (gen exacto,
    "none", o vacio/NA si no hay referencia)."""
    value = str(resultado_conocido).strip()
    lowered = value.lower()
    if any(lowered.startswith(marker) for marker in NEGATIVE_RESULT_MARKERS):
        return "none"
    if any(lowered.startswith(marker) for marker in NO_REFERENCE_RESULT_MARKERS):
        return "NA"
    first_gene = value.split(",")[0].strip()
    if first_gene.endswith("*"):
        # Variante ambigua ("CTX-M-type ESBL*"): se usa la familia sin
        # sufijo de alelo, para que compare_to_reference.py la empareje por
        # familia contra cualquier variante real que AMRFinderPlus detecte.
        from parse_amrfinder import derive_gene_family
        return derive_gene_family(first_gene.rstrip("*").strip())
    return first_gene


def find_genome_fasta(accession: str) -> Path:
    candidates = list((GENOMES_DIR / accession).glob("*_genomic.fna"))
    if not candidates:
        raise FileNotFoundError(f"No se encontro el archivo genomico para {accession} en {GENOMES_DIR / accession}")
    return candidates[0]


def run_and_log(sample_id: str, command: list[str], log_path: Path, env: dict[str, str] | None = None) -> None:
    with open(log_path, "a") as log_handle:
        log_handle.write(f"$ {shlex.join(command)}\n")
        completed = subprocess.run(command, cwd=REPO_ROOT, stdout=log_handle, stderr=subprocess.STDOUT, env=env)
    if completed.returncode != 0:
        raise RuntimeError(f"'{shlex.join(command)}' fallo con codigo {completed.returncode} (ver {log_path})")


def run_to_file(sample_id: str, command: list[str], destination: Path, log_path: Path, env: dict[str, str] | None = None) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a") as log_handle, open(destination, "w") as out_handle:
        log_handle.write(f"$ {shlex.join(command)} > {destination}\n")
        completed = subprocess.run(command, cwd=REPO_ROOT, stdout=out_handle, stderr=log_handle, env=env)
    if completed.returncode != 0:
        raise RuntimeError(f"'{shlex.join(command)}' fallo con codigo {completed.returncode} (ver {log_path})")


def process_sample(sample_id: str, accession: str, expected_genes: str, threads: int, config: dict) -> tuple[str, str]:
    tables = RESULTS_DIR / "tables"
    final_table = tables / "amr_classified" / f"{sample_id}.tsv"
    if final_table.is_file():
        return sample_id, "ya procesada, se salta"

    log_path = RESULTS_DIR / "logs" / f"{sample_id}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("")  # limpiar log de intentos previos incompletos

    try:
        genome_path = find_genome_fasta(accession)

        quast_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/quast.yaml"))
        checkm_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/checkm.yaml"))
        amrfinder_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/amrfinder.yaml"))
        abricate_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/abricate.yaml"))
        mlst_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/mlst.yaml"))

        # --- QUAST ---
        quast_dir = RESULTS_DIR / "quast" / sample_id
        quast_dir.mkdir(parents=True, exist_ok=True)
        run_and_log(sample_id, [
            "quast.py", str(genome_path), "--output-dir", str(quast_dir), "--threads", str(threads),
        ], log_path, env=quast_env)
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_quast.py"), "parse", sample_id,
            str(quast_dir / "report.tsv"), "--output-dir", str(tables / "quast"),
            "--maximum-contigs", str(config["assembly"]["maximum_contigs"]),
            "--minimum-total-length", str(config["assembly"]["minimum_total_length"]),
            "--maximum-total-length", str(config["assembly"]["maximum_total_length"]),
            "--n50-warning-threshold", str(config["assembly"]["n50_warning_threshold"]),
        ], log_path)

        # --- CheckM ---
        checkm_bin_dir = RESULTS_DIR / "checkm" / sample_id / "bins"
        checkm_out_dir = RESULTS_DIR / "checkm" / sample_id / "output"
        checkm_bin_dir.mkdir(parents=True, exist_ok=True)
        checkm_report = RESULTS_DIR / "checkm" / sample_id / "checkm_summary.tsv"
        import shutil
        shutil.copy(genome_path, checkm_bin_dir / f"{sample_id}.fasta")
        run_and_log(sample_id, ["checkm", "data", "setRoot", config["paths"]["checkm_database"]], log_path, env=checkm_env)
        run_and_log(sample_id, [
            "checkm", "lineage_wf", "-x", "fasta", "--tab_table", "-f", str(checkm_report),
            "-t", str(threads), str(checkm_bin_dir), str(checkm_out_dir),
        ], log_path, env=checkm_env)
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_checkm.py"), "parse", sample_id, str(checkm_report),
            "--output-dir", str(tables / "checkm"),
            "--minimum-completeness", str(config["assembly"]["minimum_completeness"]),
            "--maximum-contamination", str(config["assembly"]["maximum_contamination"]),
        ], log_path)

        # --- AMRFinderPlus ---
        amrfinder_table = RESULTS_DIR / "amr" / f"{sample_id}.tsv"
        amrfinder_table.parent.mkdir(parents=True, exist_ok=True)
        run_and_log(sample_id, [
            "amrfinder", "--nucleotide", str(genome_path), "--organism", "Escherichia",
            "--threads", str(threads), "--output", str(amrfinder_table),
        ], log_path, env=amrfinder_env)
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_amrfinder.py"), "parse", sample_id, str(amrfinder_table),
            "--output-dir", str(tables / "amr"),
            "--minimum-identity", str(config["amr"]["minimum_identity"]),
            "--minimum-gene-coverage", str(config["amr"]["minimum_gene_coverage"]),
        ], log_path)
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/classify_cephalosporin_genes.py"),
            str(tables / "amr" / f"{sample_id}.tsv"),
            "--resistance-targets", str(REPO_ROOT / config["resistance_targets"]),
            "--output", str(final_table),
        ], log_path)

        # --- ABricate ---
        abricate_raw_paths = []
        for database in config["amr"]["abricate_databases"]:
            raw_path = RESULTS_DIR / "abricate" / f"{sample_id}_{database}.tsv"
            run_to_file(sample_id, [
                "abricate", "--db", database, "--threads", str(threads), str(genome_path),
            ], raw_path, log_path, env=abricate_env)
            abricate_raw_paths.append(str(raw_path))
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_abricate.py"), "parse", sample_id, *abricate_raw_paths,
            "--output-dir", str(tables / "abricate"),
            "--minimum-identity", str(config["amr"]["minimum_identity"]),
            "--minimum-gene-coverage", str(config["amr"]["minimum_gene_coverage"]),
        ], log_path)

        # --- MLST ---
        mlst_table = RESULTS_DIR / "mlst" / f"{sample_id}.tsv"
        run_to_file(sample_id, [
            "mlst", "--scheme", config["mlst"]["scheme"], "--threads", str(threads), str(genome_path),
        ], mlst_table, log_path, env=mlst_env)
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/parse_mlst.py"), "parse", sample_id, str(mlst_table),
            "--output-dir", str(tables / "mlst"),
        ], log_path)

        # --- comparacion con el estandar de referencia (expected_genes) ---
        samples_row_path = RESULTS_DIR / "samples_rows" / f"{sample_id}.tsv"
        samples_row_path.parent.mkdir(parents=True, exist_ok=True)
        samples_row_path.write_text(
            "sample_id\texpected_genes\n" + f"{sample_id}\t{expected_genes}\n"
        )
        run_and_log(sample_id, [
            sys.executable, str(REPO_ROOT / "workflow/scripts/compare_to_reference.py"),
            "--samples", str(samples_row_path),
            "--amr-table", str(tables / "amr" / f"{sample_id}.tsv"),
            "--output", str(tables / "reference_comparison" / f"{sample_id}.tsv"),
        ], log_path)

        return sample_id, "OK"
    except Exception as error:  # noqa: BLE001 -- se registra en el log de la muestra, no se deja morir el hilo en silencio
        with open(log_path, "a") as log_handle:
            log_handle.write(f"ERROR: {error}\n")
        return sample_id, f"FALLO: {error}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Correr el pipeline (modo solo-ensamblaje) por lotes sobre las muestras de validacion.")
    parser.add_argument("--workers", type=int, default=4, help="Muestras procesadas en paralelo")
    parser.add_argument("--threads-per-sample", type=int, default=2, help="Hilos por herramienta, por muestra")
    parser.add_argument("--limit", type=int, default=None, help="Procesar solo las primeras N muestras (para pruebas)")
    args = parser.parse_args()

    import yaml
    with open(REPO_ROOT / "config" / "config.yaml") as config_file:
        config = yaml.safe_load(config_file)

    samples_df = pd.read_csv(SAMPLES_CSV)
    if args.limit:
        samples_df = samples_df.head(args.limit)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    jobs = []
    for _, row in samples_df.iterrows():
        accession = row["Assembly Accession"]
        sample_id = accession.replace(".", "_")
        expected_genes = derive_expected_genes(row["Resultado conocido"])
        jobs.append((sample_id, accession, expected_genes))

    print(f"Procesando {len(jobs)} muestras con {args.workers} en paralelo...", flush=True)
    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_sample, sample_id, accession, expected_genes, args.threads_per_sample, config): sample_id
            for sample_id, accession, expected_genes in jobs
        }
        for future in as_completed(futures):
            sample_id, status = future.result()
            completed_count += 1
            print(f"[{completed_count}/{len(jobs)}] {sample_id}: {status}", flush=True)


if __name__ == "__main__":
    main()
