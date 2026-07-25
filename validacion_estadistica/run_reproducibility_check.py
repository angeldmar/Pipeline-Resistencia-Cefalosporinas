"""Coeficiente de variacion (reproducibilidad) de 5 muestras de validacion.

muestras_validacion_completo.csv marca 5 filas con "control de
reproducibilidad: correr en 3 corridas independientes" en Observaciones,
pero run_validation_batch.py corre cada muestra una sola vez y no mide
tiempo/RAM -- esa parte del plan original habia quedado pendiente. Este
script la completa: corre el pipeline (modo solo-ensamblaje, mismos pasos y
mismo criterio que process_sample() en run_validation_batch.py) 3 veces
para cada una de esas 5 muestras, midiendo tiempo transcurrido y RAM maxima
por corrida (igual metodologia que workflow/scripts/run_with_timing.py:
resource.getrusage(RUSAGE_CHILDREN), portable Linux/macOS).

Escribe en resultados/reproducibilidad/{sample_id}/run{1,2,3}/, una copia
independiente de resultados/ -- nunca toca ni sobreescribe la corrida
"oficial" unica de resultados/tables/ que ya alimenta validation_summary.tsv,
para no arriesgar los resultados de la validacion principal ya confirmados.

Se corre secuencial (no en paralelo): resource.getrusage(RUSAGE_CHILDREN)
acumula el uso de TODOS los procesos hijos del proceso Python actual desde
que arranco, asi que correr replicas en paralelo dentro del mismo proceso
mezclaria sus mediciones. 15 corridas totales (5 muestras x 3), a un ritmo
similar al de run_validation_batch.py, es un tiempo razonable en secuencial.

Reanudable: si ya existe la fila de una muestra+corrida en
reproducibility_runs.tsv, esa corrida se salta.

Uso:
    python run_reproducibility_check.py [--threads N]
Salida:
    resultados/tables/reproducibility_runs.tsv (sample_id, run_number,
    elapsed_seconds, max_ram_gb -- lo que espera run_statistics.R para CV)
"""

from __future__ import annotations

import argparse
import csv
import platform
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_validation_batch import (  # noqa: E402 -- import tras sys.path.insert, patron ya usado en el resto del proyecto
    REPO_ROOT,
    SAMPLES_CSV,
    derive_expected_genes,
    env_with_conda_bin,
    find_genome_fasta,
    resolve_conda_env_bin,
    run_and_log,
    run_to_file,
)

VALIDATION_DIR = Path(__file__).resolve().parent
REPRODUCIBILITY_DIR = VALIDATION_DIR / "resultados" / "reproducibilidad"
RUNS_LOG_PATH = VALIDATION_DIR / "resultados" / "tables" / "reproducibility_runs.tsv"
RUNS_LOG_COLUMNS = ["sample_id", "run_number", "elapsed_seconds", "max_ram_gb"]

REPRODUCIBILITY_MARKER = "control de reproducibilidad"
NUMBER_OF_RUNS = 3

# Identico ajuste que run_with_timing.py: ru_maxrss se reporta en KB en
# Linux pero en BYTES en macOS/BSD.
RU_MAXRSS_IS_IN_BYTES_ON_THIS_PLATFORM = platform.system() == "Darwin"


def find_reproducibility_samples() -> list[tuple[str, str, str]]:
    """Devuelve (sample_id, accession, expected_genes) para las filas
    marcadas como control de reproducibilidad en el CSV de muestras."""
    samples_df = pd.read_csv(SAMPLES_CSV)
    flagged = samples_df[samples_df["Observaciones"].astype(str).str.contains(REPRODUCIBILITY_MARKER, na=False)]
    result = []
    for _, row in flagged.iterrows():
        accession = row["Assembly Accession"]
        sample_id = accession.replace(".", "_")
        expected_genes = derive_expected_genes(row["Resultado conocido"])
        result.append((sample_id, accession, expected_genes))
    return result


def already_logged(sample_id: str, run_number: int) -> bool:
    if not RUNS_LOG_PATH.is_file():
        return False
    existing = pd.read_csv(RUNS_LOG_PATH, sep="\t")
    return ((existing["sample_id"] == sample_id) & (existing["run_number"] == run_number)).any()


def append_run_row(sample_id: str, run_number: int, elapsed_seconds: float, max_ram_gb: float) -> None:
    RUNS_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    file_exists = RUNS_LOG_PATH.is_file()
    with open(RUNS_LOG_PATH, "a", newline="") as log_file:
        writer = csv.DictWriter(log_file, fieldnames=RUNS_LOG_COLUMNS, delimiter="\t")
        if not file_exists:
            writer.writeheader()
        writer.writerow({
            "sample_id": sample_id, "run_number": run_number,
            "elapsed_seconds": elapsed_seconds, "max_ram_gb": max_ram_gb,
        })


def run_pipeline_once(
    sample_id: str, accession: str, expected_genes: str, run_dir: Path, threads: int, config: dict,
) -> None:
    """Reproduce los mismos pasos que process_sample() de
    run_validation_batch.py (QUAST, CheckM, AMRFinderPlus, ABricate, MLST,
    comparacion con el estandar de referencia), pero escribiendo bajo
    run_dir en vez del RESULTS_DIR global de ese script -- asi cada corrida
    de reproducibilidad queda aislada y no pisa la corrida oficial unica."""
    tables = run_dir / "tables"
    log_path = run_dir / "run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("")

    genome_path = find_genome_fasta(accession)

    quast_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/quast.yaml"))
    checkm_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/checkm.yaml"))
    amrfinder_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/amrfinder.yaml"))
    abricate_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/abricate.yaml"))
    mlst_env = env_with_conda_bin(resolve_conda_env_bin("workflow/envs/mlst.yaml"))

    # --- QUAST ---
    quast_dir = run_dir / "quast" / sample_id
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
    checkm_bin_dir = run_dir / "checkm" / sample_id / "bins"
    checkm_out_dir = run_dir / "checkm" / sample_id / "output"
    checkm_bin_dir.mkdir(parents=True, exist_ok=True)
    checkm_report = run_dir / "checkm" / sample_id / "checkm_summary.tsv"
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
    amrfinder_table = run_dir / "amr" / f"{sample_id}.tsv"
    amrfinder_table.parent.mkdir(parents=True, exist_ok=True)
    final_table = tables / "amr_classified" / f"{sample_id}.tsv"
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
        raw_path = run_dir / "abricate" / f"{sample_id}_{database}.tsv"
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
    mlst_table = run_dir / "mlst" / f"{sample_id}.tsv"
    run_to_file(sample_id, [
        "mlst", "--scheme", config["mlst"]["scheme"], "--threads", str(threads), str(genome_path),
    ], mlst_table, log_path, env=mlst_env)
    run_and_log(sample_id, [
        sys.executable, str(REPO_ROOT / "workflow/scripts/parse_mlst.py"), "parse", sample_id, str(mlst_table),
        "--output-dir", str(tables / "mlst"),
    ], log_path)

    # --- comparacion con el estandar de referencia ---
    samples_row_path = run_dir / "samples_row.tsv"
    samples_row_path.write_text("sample_id\texpected_genes\n" + f"{sample_id}\t{expected_genes}\n")
    run_and_log(sample_id, [
        sys.executable, str(REPO_ROOT / "workflow/scripts/compare_to_reference.py"),
        "--samples", str(samples_row_path),
        "--amr-table", str(tables / "amr" / f"{sample_id}.tsv"),
        "--output", str(tables / "reference_comparison" / f"{sample_id}.tsv"),
    ], log_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Correr 3 veces cada muestra marcada como control de reproducibilidad, midiendo tiempo/RAM."
    )
    parser.add_argument("--threads", type=int, default=2, help="Hilos por herramienta")
    parser.add_argument(
        "--run-single", nargs=3, metavar=("SAMPLE_ID", "ACCESSION", "RUN_NUMBER"),
        help="Uso interno: corre UNA sola replica en este proceso (asi lo invoca el propio script como subproceso)",
    )
    args = parser.parse_args()

    if args.run_single:
        # Este proceso arranca fresco (lo lanza el bloque de abajo via
        # subprocess.run, uno nuevo por replica), asi que
        # resource.getrusage(RUSAGE_CHILDREN) parte de cero: una sola
        # lectura al final, SIN restar un "antes", ya refleja exactamente el
        # pico de RAM de sus hijos (quast/checkm/amrfinder/abricate/mlst) en
        # ESTA replica. Restar dos lecturas tomadas del mismo proceso de
        # larga vida (como se intento primero, midiendo desde el bucle de
        # abajo) NO funciona: ru_maxrss es un maximo acumulado que nunca
        # baja, asi que una replica con picos mas bajos que una anterior
        # quedaria subestimada o en cero. Mismo principio que ya usa
        # run_with_timing.py (una invocacion == un proceso == una medicion limpia).
        sample_id, accession, run_number_str = args.run_single
        run_number = int(run_number_str)
        with open(REPO_ROOT / "config" / "config.yaml") as config_file:
            config = yaml.safe_load(config_file)
        samples_df = pd.read_csv(SAMPLES_CSV)
        expected_genes = derive_expected_genes(
            samples_df.loc[samples_df["Assembly Accession"] == accession, "Resultado conocido"].iloc[0]
        )
        run_dir = REPRODUCIBILITY_DIR / sample_id / f"run{run_number}"
        run_dir.mkdir(parents=True, exist_ok=True)

        start_wall_clock = time.perf_counter()
        run_pipeline_once(sample_id, accession, expected_genes, run_dir, args.threads, config)
        elapsed_seconds = round(time.perf_counter() - start_wall_clock, 3)

        resource_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        max_rss_kilobytes = resource_usage.ru_maxrss / 1024 if RU_MAXRSS_IS_IN_BYTES_ON_THIS_PLATFORM else resource_usage.ru_maxrss
        max_ram_gb = round(max_rss_kilobytes / (1024 * 1024), 4)

        append_run_row(sample_id, run_number, elapsed_seconds, max_ram_gb)
        print(f"{sample_id} corrida {run_number}: OK ({elapsed_seconds}s, {max_ram_gb} GB)", flush=True)
        return

    samples = find_reproducibility_samples()
    print(f"{len(samples)} muestra(s) de control de reproducibilidad, {NUMBER_OF_RUNS} corridas cada una.", flush=True)

    for sample_id, accession, _expected_genes in samples:
        for run_number in range(1, NUMBER_OF_RUNS + 1):
            if already_logged(sample_id, run_number):
                print(f"{sample_id} corrida {run_number}: ya registrada, se salta", flush=True)
                continue

            # Cada replica se lanza como un PROCESO NUEVO -- ver comentario
            # en la rama --run-single de arriba sobre por que la medicion
            # de RAM tiene que ocurrir alli adentro, no aqui.
            log_path = REPRODUCIBILITY_DIR / sample_id / f"run{run_number}.log"
            log_path.parent.mkdir(parents=True, exist_ok=True)

            with open(log_path, "w") as log_file:
                completed = subprocess.run(
                    [sys.executable, __file__, "--run-single", sample_id, accession, str(run_number),
                     "--threads", str(args.threads)],
                    cwd=REPO_ROOT, stdout=log_file, stderr=subprocess.STDOUT,
                )

            if completed.returncode != 0:
                print(f"{sample_id} corrida {run_number}: FALLO (ver {log_path})", flush=True)
            else:
                # El propio subproceso --run-single ya escribio su fila en
                # reproducibility_runs.tsv y su resumen a log_path.
                with open(log_path) as log_file:
                    last_line = log_file.readlines()[-1].strip() if log_file else ""
                print(last_line or f"{sample_id} corrida {run_number}: OK", flush=True)


if __name__ == "__main__":
    main()
