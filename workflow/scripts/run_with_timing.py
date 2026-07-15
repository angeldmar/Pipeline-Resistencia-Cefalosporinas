"""Envoltorio (wrapper) que mide tiempo y memoria de un comando externo.

El diseno del pipeline sugiere usar "/usr/bin/time -v" para medir tiempo y
RAM de cada herramienta. Esa variante (-v, "verbose") es especifica de GNU
time y no existe en macOS/BSD ni se puede asumir instalada en cualquier
entorno de ejecucion. En su lugar, este script logra lo mismo con el modulo
"resource" de la biblioteca estandar de Python (portable entre Linux y
macOS), manteniendo el principio ya usado en todo el pipeline de que Python
orquesta las herramientas externas via subprocess, en vez de depender de
utilidades de shell especificas de una plataforma.

Uso: se antepone a cualquier comando real, separado por "--":

    python workflow/scripts/run_with_timing.py \
      --sample-id EC001 --module spades --threads 8 \
      --output results/tables/performance/EC001_spades.tsv \
      -- spades.py -1 R1.fastq.gz -2 R2.fastq.gz -o outdir -t 8 --careful

El comando real se ejecuta tal cual (su stdout/stderr se heredan del
wrapper, asi que la redireccion ">{log} 2>&1" que ya usan las reglas de
Snakemake lo sigue capturando sin cambios). El codigo de salida del comando
real se propaga como codigo de salida de este script, para que Snakemake
detecte correctamente si la regla fallo.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import platform
import resource
import subprocess
import sys
import time

import pandas as pd

PERFORMANCE_LOG_COLUMNS = [
    "sample_id",
    "module",
    "elapsed_seconds",
    "cpu_seconds",
    "max_ram_gb",
    "exit_code",
    "threads",
    "run_date",
]

# resource.getrusage().ru_maxrss se reporta en KILOBYTES en Linux, pero en
# BYTES en macOS/BSD. Sin este ajuste, la RAM maxima quedaria sobrestimada
# por un factor de ~1000 en macOS.
RU_MAXRSS_IS_IN_BYTES_ON_THIS_PLATFORM = platform.system() == "Darwin"


def run_command_with_timing(command: list[str]) -> tuple[int, float, float, float]:
    """Ejecuta command y devuelve (codigo_de_salida, segundos_transcurridos,
    segundos_de_cpu, RAM_maxima_en_GB) del proceso hijo."""
    start_wall_clock = time.perf_counter()
    completed_process = subprocess.run(command, check=False)
    elapsed_seconds = time.perf_counter() - start_wall_clock

    # RUSAGE_CHILDREN acumula el uso de recursos de los procesos hijos
    # terminados hasta ahora; como este script es un proceso nuevo por cada
    # invocacion (una por regla de Snakemake), esto refleja unicamente al
    # comando que se acaba de ejecutar.
    resource_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_seconds = round(resource_usage.ru_utime + resource_usage.ru_stime, 3)

    max_rss_kilobytes = (
        resource_usage.ru_maxrss / 1024 if RU_MAXRSS_IS_IN_BYTES_ON_THIS_PLATFORM else resource_usage.ru_maxrss
    )
    max_ram_gb = round(max_rss_kilobytes / (1024 * 1024), 4)

    return completed_process.returncode, round(elapsed_seconds, 3), cpu_seconds, max_ram_gb


def write_performance_row(
    output_path: Path,
    sample_id: str,
    module: str,
    threads: int,
    exit_code: int,
    elapsed_seconds: float,
    cpu_seconds: float,
    max_ram_gb: float,
) -> None:
    """Escribe una fila con las metricas de desempeno de esta ejecucion."""
    performance_row = {
        "sample_id": sample_id,
        "module": module,
        "elapsed_seconds": elapsed_seconds,
        "cpu_seconds": cpu_seconds,
        "max_ram_gb": max_ram_gb,
        "exit_code": exit_code,
        "threads": threads,
        "run_date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([performance_row], columns=PERFORMANCE_LOG_COLUMNS).to_csv(output_path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ejecutar un comando midiendo tiempo, CPU y RAM, y registrar el resultado.",
    )
    parser.add_argument("--sample-id", required=True, help="Identificador de la muestra que se esta procesando")
    parser.add_argument("--module", required=True, help="Nombre del modulo/herramienta que se esta midiendo (ej. spades, fastp)")
    parser.add_argument("--threads", type=int, default=1, help="Numero de hilos con los que se invoco el comando")
    parser.add_argument("--output", type=Path, required=True, help="Ruta del TSV de desempeno de esta ejecucion")
    parser.add_argument(
        "command", nargs=argparse.REMAINDER,
        help="Comando real a ejecutar, precedido por '--' (ej. -- spades.py -1 ... )",
    )
    args = parser.parse_args()

    # argparse.REMAINDER conserva el "--" separador si el usuario lo puso;
    # se descarta para quedarnos solo con el comando real.
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("Falta el comando a ejecutar despues de '--'")

    exit_code, elapsed_seconds, cpu_seconds, max_ram_gb = run_command_with_timing(command)

    write_performance_row(
        args.output, args.sample_id, args.module, args.threads,
        exit_code, elapsed_seconds, cpu_seconds, max_ram_gb,
    )

    print(
        f"[{args.module}/{args.sample_id}] tiempo: {elapsed_seconds}s, CPU: {cpu_seconds}s, "
        f"RAM maxima: {max_ram_gb} GB, codigo de salida: {exit_code}"
    )

    # Propaga el codigo de salida del comando real, para que Snakemake
    # detecte correctamente si la regla fallo (nunca se debe "tragar" un
    # error solo porque ya se registro la metrica de desempeno).
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
