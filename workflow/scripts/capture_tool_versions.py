"""Registra las versiones instaladas de las herramientas externas del pipeline.

A diferencia de las metricas de desempeno (una fila por muestra y modulo),
la version de una herramienta no cambia entre muestras dentro de una misma
corrida del pipeline: es una propiedad del entorno, no de los datos. Por eso
este script se corre UNA sola vez por corrida (no por muestra) y deja un
unico archivo de trazabilidad: data/metadata/tool_versions.tsv.

Cada reporte individual (generate_report.py) lee este archivo para mostrar
con que version de cada herramienta se generaron sus resultados.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import subprocess

import pandas as pd

# Comando usado para consultar la version de cada herramienta externa del
# pipeline. Snakemake ya fija la version exacta de cada una via su ambiente
# conda (workflow/envs/*.yaml); esto solo confirma, en tiempo de ejecucion,
# cual version quedo realmente instalada.
TOOL_VERSION_COMMANDS = {
    "fastp": ["fastp", "--version"],
    "spades": ["spades.py", "--version"],
    "quast": ["quast.py", "--version"],
    "checkm": ["checkm", "--version"],
    "kraken2": ["kraken2", "--version"],
    "prokka": ["prokka", "--version"],
    "amrfinder": ["amrfinder", "--version"],
    "snakemake": ["snakemake", "--version"],
}


def get_tool_version(command: list[str]) -> str:
    """Ejecuta el comando de version de una herramienta y devuelve su
    salida. Si la herramienta no esta instalada en el entorno actual (por
    ejemplo, al correr este script fuera de sus ambientes conda), devuelve
    "not installed" en vez de fallar: la version es un dato de trazabilidad,
    no debe bloquear el resto del pipeline."""
    try:
        completed_process = subprocess.run(command, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return "not installed"

    # Distintas herramientas imprimen la version en stdout o en stderr.
    version_output = (completed_process.stdout + completed_process.stderr).strip()
    return version_output.splitlines()[0] if version_output else "unknown"


def capture_all_tool_versions() -> pd.DataFrame:
    """Consulta la version de cada herramienta del pipeline y las junta en
    una tabla, junto con la fecha en que se consultaron."""
    captured_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    version_rows = [
        {"tool": tool_name, "version": get_tool_version(command), "captured_date": captured_date}
        for tool_name, command in TOOL_VERSION_COMMANDS.items()
    ]
    return pd.DataFrame(version_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Registrar las versiones instaladas de las herramientas externas del pipeline."
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/metadata/tool_versions.tsv"),
    )
    args = parser.parse_args()

    tool_versions_table = capture_all_tool_versions()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tool_versions_table.to_csv(args.output, sep="\t", index=False)

    print(f"Versiones de {len(tool_versions_table)} herramienta(s) registradas en {args.output}")


if __name__ == "__main__":
    main()
