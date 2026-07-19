"""Orquestacion del pipeline para la interfaz web local de analisis ad-hoc.

Esta interfaz es una herramienta de conveniencia DISTINTA del pipeline en si:
el pipeline (Snakemake + workflow/) esta pensado para lotes curados de
muestras publicas, documentadas en config/samples.tsv. Esta interfaz existe
para el caso de uso "tengo un FASTQ o un FASTA suelto y quiero ver que
produce el pipeline con el", sin tener que editar samples.tsv a mano.

Dos modos de entrada, con dos formas de orquestacion distintas:

  FASTQ (lecturas crudas paired-end): modo COMPLETO. Se arma una tabla de
  muestras AISLADA de una sola fila (nunca se toca config/samples.tsv, para
  no mezclar una carga ad-hoc con el lote curado principal) y se lanza
  Snakemake de verdad (`--use-conda`), apuntando al reporte completo. Corre
  el pipeline entero: fastp, SPAdes, QUAST, CheckM, Kraken2, AMRFinderPlus,
  ABricate, MLST, Prokka.

  FASTA (ensamblaje ya armado): modo PARCIAL. No hay lecturas crudas, asi
  que fastp/cobertura/taxonomia no pueden correr (Kraken2 necesita lecturas,
  no un ensamblaje ya colapsado). En vez de forzar a Snakemake a aceptar un
  grafo de dependencias incompleto (fragil: se probo y el comportamiento de
  "--touch" para simular pasos previos no fue confiable de un intento a
  otro), este modo ejecuta secuencialmente, con subprocess, las MISMAS
  herramientas y con las MISMAS rutas de archivo que usarian las reglas de
  Snakemake (QUAST, CheckM, Prokka, AMRFinderPlus, ABricate, MLST), y
  reutiliza sin modificar merge_results.py y generate_report.py -- ambos ya
  toleran modulos faltantes (fastp/taxonomia quedan "N/D" en el reporte,
  mismo mecanismo de resiliencia ya probado desde la parte 17).
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import os
import re
import shlex
import subprocess
import sys
import threading

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = Path(__file__).resolve().parent / "runs"
SNAKEMAKE_CONDA_DIR = REPO_ROOT / ".snakemake" / "conda"

# Solo letras, numeros, guiones y guiones bajos: evita traversal de rutas
# (../) y caracteres que romperian nombres de archivo o comandos de shell.
SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

# Estas muestras nunca vienen de un repositorio publico documentado, asi que
# no tienen una accesion SRA/BioSample real. Se usa un marcador fijo, valido
# en formato (para no romper validate_samples.py) pero inconfundible: el
# campo "data_source" es el que realmente deja constancia de que la muestra
# es una carga local ad-hoc, no una accesion real.
PLACEHOLDER_RUN_ACCESSION = "SRR000000"
PLACEHOLDER_BIOSAMPLE = "SAMN00000000"
WEBAPP_DATA_SOURCE = "Carga local ad-hoc (interfaz web)"


class MissingCondaEnvironmentError(RuntimeError):
    pass


def resolve_conda_env_bin(env_yaml_relative_path: str) -> Path:
    """Ubica el bin/ del ambiente Conda que Snakemake ya creo para la regla
    correspondiente (workflow/envs/<tool>.yaml), sin reimplementar el
    algoritmo de hash interno de Snakemake (no es una API publica ni
    estable entre versiones). En vez de eso, compara el CONTENIDO del yaml
    de origen contra los archivos marcador que Snakemake deja en
    .snakemake/conda/<hash>.yaml -- que son una copia textual exacta del
    yaml que los genero -- y devuelve el ambiente cuyo .env_setup_done
    confirma que la instalacion termino, no solo que empezo."""
    source_content = (REPO_ROOT / env_yaml_relative_path).read_text()

    if SNAKEMAKE_CONDA_DIR.is_dir():
        for marker in SNAKEMAKE_CONDA_DIR.glob("*.yaml"):
            if marker.read_text() == source_content:
                env_hash = marker.stem
                bin_dir = SNAKEMAKE_CONDA_DIR / env_hash / "bin"
                if (SNAKEMAKE_CONDA_DIR / f"{env_hash}.env_setup_done").is_file() and bin_dir.is_dir():
                    return bin_dir

    raise MissingCondaEnvironmentError(
        f"No se encontro un ambiente Conda ya creado para {env_yaml_relative_path}. "
        "Crea los ambientes del pipeline primero, por ejemplo con: "
        "CONDA_SUBDIR=osx-64 snakemake --use-conda --conda-create-envs-only --cores 1"
    )


def _env_with_conda_bin(bin_dir: Path) -> dict[str, str]:
    """Variables de entorno para invocar una herramienta de un ambiente
    Conda de Snakemake por subprocess directo (sin pasar por Snakemake ni
    'conda activate'): antepone el bin/ del ambiente al PATH actual, para
    que tanto el ejecutable principal como sus dependencias empaquetadas en
    el mismo ambiente (ej. blastn para abricate, prodigal/hmmer para
    checkm) se encuentren, igual que si el ambiente estuviera activado."""
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    return env


class InvalidSampleIdError(ValueError):
    pass


def validate_sample_id(sample_id: str) -> None:
    if not SAMPLE_ID_PATTERN.match(sample_id):
        raise InvalidSampleIdError(
            "El identificador de muestra solo puede contener letras, numeros, "
            "guiones y guiones bajos (maximo 64 caracteres)."
        )


def job_dir(sample_id: str) -> Path:
    """Carpeta de trabajo de esta carga: samples.tsv aislado, log y estado."""
    directory = RUNS_DIR / sample_id
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def status_file(sample_id: str) -> Path:
    return job_dir(sample_id) / "status.txt"


def log_file(sample_id: str) -> Path:
    return job_dir(sample_id) / "run.log"


def write_status(sample_id: str, status: str) -> None:
    status_file(sample_id).write_text(status)


def read_status(sample_id: str) -> str:
    path = status_file(sample_id)
    return path.read_text().strip() if path.is_file() else "not_found"


def report_path(sample_id: str) -> Path:
    return REPO_ROOT / "results" / "reports" / f"{sample_id}.html"


def create_isolated_samples_file(sample_id: str, sequencing_platform: str, expected_genes: str) -> Path:
    """Arma una tabla de muestras de una sola fila para esta carga, separada
    de config/samples.tsv. expected_genes queda "NA" si no se especifica: el
    pipeline ya interpreta eso como "sin estandar de referencia documentado"
    (categoria "Indeterminado" en la comparacion, ver compare_to_reference.py),
    en vez de forzar una comparacion sin sentido contra un valor inventado."""
    samples_path = job_dir(sample_id) / "samples.tsv"
    header = [
        "sample_id", "run_accession", "biosample", "sequencing_platform",
        "phenotype_cefotaxime", "phenotype_ceftriaxone", "phenotype_ceftazidime",
        "expected_genes", "data_source",
    ]
    row = [
        sample_id, PLACEHOLDER_RUN_ACCESSION, PLACEHOLDER_BIOSAMPLE, sequencing_platform or "ILLUMINA",
        "NA", "NA", "NA", expected_genes.strip() or "NA", WEBAPP_DATA_SOURCE,
    ]
    samples_path.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n")
    return samples_path


def save_uploaded_fastq(sample_id: str, r1_storage, r2_storage) -> None:
    """Guarda el par de FASTQ subidos en data/raw/, con el nombre exacto que
    el pipeline espera (ver download_sample / fastp)."""
    raw_dir = REPO_ROOT / "data" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    r1_storage.save(raw_dir / f"{sample_id}_R1.fastq.gz")
    r2_storage.save(raw_dir / f"{sample_id}_R2.fastq.gz")


def save_uploaded_fasta(sample_id: str, fasta_storage) -> Path:
    """Guarda el FASTA subido directamente como el ensamblaje ya filtrado,
    en la misma ruta que usaria filter_contigs.py (results/assemblies/.../
    contigs.filtered.fasta), para que las herramientas de QC/AMR/tipificacion
    lo consuman sin ningun paso intermedio."""
    assembly_dir = REPO_ROOT / "results" / "assemblies" / sample_id
    assembly_dir.mkdir(parents=True, exist_ok=True)
    assembly_path = assembly_dir / "contigs.filtered.fasta"
    fasta_storage.save(assembly_path)
    return assembly_path


def _touch_raw_fastq_files(sample_id: str) -> None:
    """Actualiza la fecha de modificacion de los FASTQ ya subidos a "ahora".
    create_isolated_samples_file() se llama DESPUES de guardar esos FASTQ
    (ver launch_fastq_pipeline), asi que sin esto, samples.tsv (entrada de
    la regla download_sample) queda con una fecha mas reciente que sus
    propias salidas (data/raw/{sample}_R1/R2.fastq.gz) -- Snakemake
    interpretaria eso como "salida desactualizada" y reintentaria descargar
    la accesion SRA placeholder (SRR000000, invalida) en vez de usar los
    FASTQ recien subidos."""
    raw_dir = REPO_ROOT / "data" / "raw"
    for suffix in ("R1", "R2"):
        path = raw_dir / f"{sample_id}_{suffix}.fastq.gz"
        if path.is_file():
            path.touch()


def _write_placeholder_download_performance(sample_id: str) -> None:
    """Crea la tercera salida declarada de la regla download_sample
    (el registro de desempeno de la "descarga", en ceros) para que
    Snakemake la vea como ya satisfecha y NO reintente descargar la
    accesion SRA placeholder. Sin esto, aunque los FASTQ ya esten en su
    lugar (ver _touch_raw_fastq_files), Snakemake igual dispara
    download_sample por faltarle esta salida -- y si esa descarga falla
    (la accesion placeholder es invalida), el mecanismo de limpieza de
    Snakemake borra TODAS las salidas declaradas de la regla fallida,
    incluidos los FASTQ reales que nunca se tocaron. Mismo esquema de
    columnas que escribe run_with_timing.py (ver PERFORMANCE_LOG_COLUMNS),
    para que combine_performance.py lo lea sin problema."""
    performance_path = REPO_ROOT / "results" / "tables" / "performance" / f"{sample_id}_download.tsv"
    performance_path.parent.mkdir(parents=True, exist_ok=True)
    columns = ["sample_id", "module", "elapsed_seconds", "cpu_seconds", "max_ram_gb", "exit_code", "threads", "run_date"]
    row = [
        sample_id, "download", "0", "0", "0", "0", "0",
        datetime.now(timezone.utc).strftime("%Y-%m-%d"),
    ]
    performance_path.write_text("\t".join(columns) + "\n" + "\t".join(row) + "\n")


def _append_log(sample_id: str, message: str) -> None:
    with open(log_file(sample_id), "a") as handle:
        handle.write(f"[{datetime.now(timezone.utc).isoformat()}] {message}\n")


def launch_fastq_pipeline(sample_id: str, sequencing_platform: str, expected_genes: str, threads: int) -> None:
    """Lanza el pipeline COMPLETO via Snakemake, en un proceso en segundo
    plano. No bloquea la peticion HTTP que la dispara."""
    samples_path = create_isolated_samples_file(sample_id, sequencing_platform, expected_genes)
    _touch_raw_fastq_files(sample_id)
    _write_placeholder_download_performance(sample_id)
    write_status(sample_id, "running")

    command = [
        "snakemake",
        "--use-conda",
        "--cores", str(threads),
        "--rerun-incomplete",
        str(report_path(sample_id).relative_to(REPO_ROOT)),
        # "--config" espera uno o mas pares clave=valor (nargs="+") y
        # consume todo lo que le siga en la linea de comandos: si el
        # archivo objetivo fuera despues, Snakemake lo interpretaria como
        # otro par clave=valor invalido. Por eso va al final.
        "--config", f"samples={samples_path}",
    ]
    _append_log(sample_id, f"Comando: {shlex.join(command)}")

    def _run() -> None:
        with open(log_file(sample_id), "a") as log_handle:
            completed = subprocess.run(command, cwd=REPO_ROOT, stdout=log_handle, stderr=subprocess.STDOUT)
        write_status(sample_id, "done" if completed.returncode == 0 else f"failed:codigo {completed.returncode}")

    threading.Thread(target=_run, daemon=True).start()


def _run_and_log(sample_id: str, command: list[str], env: dict[str, str] | None = None) -> None:
    """Corre un comando y agrega su salida al log de la carga. Lanza
    RuntimeError con un mensaje claro si el comando falla (incluida la
    herramienta no encontrada), para que el hilo de FASTA pueda detenerse
    y marcar la carga como fallida en vez de continuar con datos a medias."""
    _append_log(sample_id, f"Ejecutando: {shlex.join(command)}")
    try:
        completed = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True, env=env)
    except FileNotFoundError as error:
        raise RuntimeError(f"Herramienta no encontrada: {command[0]} ({error})") from error

    with open(log_file(sample_id), "a") as handle:
        handle.write(completed.stdout)
        handle.write(completed.stderr)

    if completed.returncode != 0:
        raise RuntimeError(f"'{shlex.join(command)}' fallo con codigo {completed.returncode}")


def run_fasta_only_pipeline(sample_id: str, assembly_path: Path, expected_genes: str, threads: int) -> None:
    """Corre, en un hilo en segundo plano, la secuencia de herramientas que
    SI pueden operar solo con un ensamblaje (sin lecturas crudas): QUAST,
    CheckM, Prokka, AMRFinderPlus, ABricate y MLST -- las mismas rutas de
    archivo que usarian las reglas de Snakemake equivalentes, para que
    merge_results.py y generate_report.py las encuentren sin cambios."""
    write_status(sample_id, "running")

    def _run() -> None:
        try:
            with open(REPO_ROOT / "config" / "config.yaml") as config_file:
                config = yaml.safe_load(config_file)

            results = REPO_ROOT / "results"
            tables = results / "tables"
            performance_dir = tables / "performance"
            for directory in (
                results / "qc" / "quast" / sample_id, results / "qc" / "checkm" / sample_id / "bins",
                results / "qc" / "checkm" / sample_id / "output", results / "annotation" / sample_id,
                tables / "quast", tables / "checkm", tables / "amr", tables / "abricate", tables / "mlst",
                performance_dir,
            ):
                directory.mkdir(parents=True, exist_ok=True)

            # Estas herramientas solo existen dentro de los ambientes Conda
            # aislados que Snakemake crea por regla (.snakemake/conda/), no
            # en el PATH general de la shell -- se resuelven aqui, una sola
            # vez, y se antepone su bin/ al PATH de cada subprocess que las
            # invoca (ver resolve_conda_env_bin / _env_with_conda_bin).
            quast_env = _env_with_conda_bin(resolve_conda_env_bin("workflow/envs/quast.yaml"))
            checkm_env = _env_with_conda_bin(resolve_conda_env_bin("workflow/envs/checkm.yaml"))
            amrfinder_env = _env_with_conda_bin(resolve_conda_env_bin("workflow/envs/amrfinder.yaml"))
            abricate_env = _env_with_conda_bin(resolve_conda_env_bin("workflow/envs/abricate.yaml"))
            mlst_env = _env_with_conda_bin(resolve_conda_env_bin("workflow/envs/mlst.yaml"))

            def timed(module: str, command: list[str], env: dict[str, str] | None = None) -> None:
                """Para herramientas con su propio flag de salida (--output,
                -f, etc.): el resultado va a un archivo declarado, y stdout/
                stderr combinados (incluido el mensaje de run_with_timing.py)
                se registran en el log de la carga sin problema."""
                _run_and_log(sample_id, [
                    sys.executable, str(REPO_ROOT / "workflow/scripts/run_with_timing.py"),
                    "--sample-id", sample_id, "--module", module, "--threads", str(threads),
                    "--output", str(performance_dir / f"{sample_id}_{module}.tsv"),
                    "--", *command,
                ], env=env)

            def timed_to_file(module: str, command: list[str], destination: Path, env: dict[str, str] | None = None) -> None:
                """Para herramientas que escriben su resultado por stdout
                (abricate, mlst): stdout debe ir EXCLUSIVAMENTE al archivo de
                destino, nunca mezclarse con el log (run_with_timing.py ya
                envia su propio mensaje de estado a stderr especificamente
                para permitir este uso, ver parte 25)."""
                _append_log(sample_id, f"Ejecutando ({module}): {shlex.join(command)} > {destination}")
                destination.parent.mkdir(parents=True, exist_ok=True)
                with open(destination, "w") as output_handle, open(log_file(sample_id), "a") as log_handle:
                    completed = subprocess.run(
                        [
                            sys.executable, str(REPO_ROOT / "workflow/scripts/run_with_timing.py"),
                            "--sample-id", sample_id, "--module", module, "--threads", str(threads),
                            "--output", str(performance_dir / f"{sample_id}_{module}.tsv"),
                            "--", *command,
                        ],
                        cwd=REPO_ROOT, stdout=output_handle, stderr=log_handle, env=env,
                    )
                if completed.returncode != 0:
                    raise RuntimeError(f"'{shlex.join(command)}' fallo con codigo {completed.returncode}")

            # --- QUAST -----------------------------------------------------
            quast_report = results / "qc" / "quast" / sample_id / "report.tsv"
            timed(
                "quast", ["quast.py", str(assembly_path), "--output-dir", str(quast_report.parent), "--threads", str(threads)],
                env=quast_env,
            )
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/parse_quast.py"), "parse", sample_id, str(quast_report),
                "--output-dir", str(tables / "quast"),
                "--maximum-contigs", str(config["assembly"]["maximum_contigs"]),
                "--minimum-total-length", str(config["assembly"]["minimum_total_length"]),
                "--maximum-total-length", str(config["assembly"]["maximum_total_length"]),
                "--n50-warning-threshold", str(config["assembly"]["n50_warning_threshold"]),
            ])

            # --- CheckM ------------------------------------------------------
            checkm_bin_dir = results / "qc" / "checkm" / sample_id / "bins"
            checkm_out_dir = results / "qc" / "checkm" / sample_id / "output"
            checkm_report = results / "qc" / "checkm" / sample_id / "checkm_summary.tsv"
            import shutil
            shutil.copy(assembly_path, checkm_bin_dir / f"{sample_id}.fasta")
            _run_and_log(sample_id, ["checkm", "data", "setRoot", config["paths"]["checkm_database"]], env=checkm_env)
            timed("checkm", [
                "checkm", "lineage_wf", "-x", "fasta", "--tab_table", "-f", str(checkm_report),
                "-t", str(threads), str(checkm_bin_dir), str(checkm_out_dir),
            ], env=checkm_env)
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/parse_checkm.py"), "parse", sample_id, str(checkm_report),
                "--output-dir", str(tables / "checkm"),
                "--minimum-completeness", str(config["assembly"]["minimum_completeness"]),
                "--maximum-contamination", str(config["assembly"]["maximum_contamination"]),
            ])

            # --- AMRFinderPlus -------------------------------------------------
            amrfinder_table = results / "amr" / "amrfinder" / f"{sample_id}.tsv"
            amrfinder_table.parent.mkdir(parents=True, exist_ok=True)
            timed("amrfinder", [
                "amrfinder", "--nucleotide", str(assembly_path), "--organism", "Escherichia",
                "--threads", str(threads), "--output", str(amrfinder_table),
            ], env=amrfinder_env)
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/parse_amrfinder.py"), "parse", sample_id, str(amrfinder_table),
                "--output-dir", str(tables / "amr"),
                "--minimum-identity", str(config["amr"]["minimum_identity"]),
                "--minimum-gene-coverage", str(config["amr"]["minimum_gene_coverage"]),
            ])
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/classify_cephalosporin_genes.py"),
                str(tables / "amr" / f"{sample_id}.tsv"),
                "--resistance-targets", str(REPO_ROOT / config["resistance_targets"]),
                "--output", str(tables / f"{sample_id}_amr_classified.tsv"),
            ])

            # --- ABricate (segundo motor) --------------------------------------
            abricate_raw_paths = []
            for database in config["amr"]["abricate_databases"]:
                raw_path = results / "amr" / "abricate" / f"{sample_id}_{database}.tsv"
                timed_to_file(
                    f"abricate_{database}",
                    ["abricate", "--db", database, "--threads", str(threads), str(assembly_path)],
                    raw_path,
                    env=abricate_env,
                )
                abricate_raw_paths.append(str(raw_path))
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/parse_abricate.py"), "parse", sample_id, *abricate_raw_paths,
                "--output-dir", str(tables / "abricate"),
                "--minimum-identity", str(config["amr"]["minimum_identity"]),
                "--minimum-gene-coverage", str(config["amr"]["minimum_gene_coverage"]),
            ])

            # --- MLST ------------------------------------------------------
            mlst_table = results / "typing" / "mlst" / f"{sample_id}.tsv"
            timed_to_file(
                "mlst",
                ["mlst", "--scheme", config["mlst"]["scheme"], "--threads", str(threads), str(assembly_path)],
                mlst_table,
                env=mlst_env,
            )
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/parse_mlst.py"), "parse", sample_id, str(mlst_table),
                "--output-dir", str(tables / "mlst"),
            ])

            # --- Concordancia entre motores de AMR (ya se tienen ambas tablas) ---
            engine_concordance_output = tables / f"{sample_id}_engine_concordance.tsv"
            samples_path = create_isolated_samples_file(sample_id, "N/A (solo ensamblaje)", expected_genes)
            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/compare_amr_engines.py"),
                "--samples", str(samples_path),
                "--amrfinder-table", str(tables / "amr" / f"{sample_id}.tsv"),
                "--abricate-table", str(tables / "abricate" / f"{sample_id}.tsv"),
                "--concordance-output", str(engine_concordance_output),
                "--agreement-output", str(job_dir(sample_id) / "engine_concordance_input.csv"),
            ])

            # --- Integracion y reporte (solo esta muestra) ----------------------
            # Los modulos que este modo NO corre (fastp, taxonomia, comparacion
            # de referencia, desempeno agregado) se apuntan a rutas que no
            # existen dentro de la carpeta de esta carga: merge_results.py los
            # omite con un aviso (ver load_optional_table) en vez de -- por
            # accidente -- recoger el archivo COMPARTIDO del pipeline principal
            # si este ya se corrio antes para las muestras curadas.
            missing_module_path = job_dir(sample_id) / "not_available.tsv"
            master_output = tables / f"{sample_id}_master_results.tsv"
            merge_command = [
                sys.executable, str(REPO_ROOT / "workflow/scripts/merge_results.py"),
                "--samples", str(samples_path),
                "--fastp-summary", str(missing_module_path),
                "--quast-summary", str(tables / "quast" / f"{sample_id}.tsv"),
                "--checkm-summary", str(tables / "checkm" / f"{sample_id}.tsv"),
                "--taxonomy-summary", str(missing_module_path),
                "--amr-summary", str(tables / "amr" / f"{sample_id}.tsv"),
                "--reference-comparison", str(missing_module_path),
                "--performance-by-sample", str(missing_module_path),
                "--engine-concordance", str(engine_concordance_output),
                "--mlst-summary", str(tables / "mlst" / f"{sample_id}.tsv"),
                "--output", str(master_output),
            ]
            _run_and_log(sample_id, merge_command)

            _run_and_log(sample_id, [
                sys.executable, str(REPO_ROOT / "workflow/scripts/generate_report.py"),
                "--sample-id", sample_id,
                "--master-results", str(master_output),
                "--amr-classified", str(tables / f"{sample_id}_amr_classified.tsv"),
                "--performance-summary", str(missing_module_path),
                "--output", str(report_path(sample_id)),
            ])

            _append_log(sample_id, "Analisis parcial (solo ensamblaje) completado.")
            write_status(sample_id, "done")
        except Exception as error:  # noqa: BLE001 -- se registra y se expone en /status, no se deja morir el hilo en silencio
            _append_log(sample_id, f"ERROR: {error}")
            write_status(sample_id, f"failed:{error}")

    threading.Thread(target=_run, daemon=True).start()
