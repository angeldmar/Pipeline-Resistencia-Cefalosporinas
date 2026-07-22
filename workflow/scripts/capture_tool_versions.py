"""Registra las versiones instaladas de las herramientas externas del pipeline.

A diferencia de las metricas de desempeno (una fila por muestra y modulo),
la version de una herramienta no cambia entre muestras dentro de una misma
corrida del pipeline: es una propiedad del entorno, no de los datos. Por eso
este script se corre UNA sola vez por corrida (no por muestra) y deja un
unico archivo de trazabilidad: data/metadata/tool_versions.tsv.

Cada reporte individual (generate_report.py) lee este archivo para mostrar
con que version de cada herramienta se generaron sus resultados.

Cada herramienta externa vive en su PROPIO ambiente Conda aislado (uno por
regla de Snakemake), nunca en el ambiente generico de este script
(workflow/envs/python.yaml) -- llamarlas por nombre suelto en el PATH
actual siempre devolveria "not installed", sin importar si de verdad estan
instaladas en su propio ambiente. Por eso cada herramienta se resuelve
explicitamente contra su propio ambiente ya creado (ver
resolve_conda_env_bin, mismo patron que webapp/pipeline_runner.py y
workflow/rules/assembly.smk). "snakemake" es la unica excepcion: no vive en
un ambiente por-regla (es quien orquesta todo), asi que se consulta
importando el paquete directamente en vez de por subprocess -- funciona
porque la regla que invoca este script usa {sys.executable} (el interprete
que lanzo Snakemake), no el "python" del ambiente activado.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import os
import re
import subprocess

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SNAKEMAKE_CONDA_DIR = REPO_ROOT / ".snakemake" / "conda"


class MissingCondaEnvironmentError(RuntimeError):
    pass


def resolve_conda_env_bin(env_yaml_relative_path: str) -> Path:
    """Ubica el bin/ del ambiente Conda ya creado por Snakemake para una
    herramienta, comparando el contenido del yaml de origen contra los
    marcadores en .snakemake/conda/ (duplicado a proposito, mismo criterio
    de scripts autocontenidos que en webapp/pipeline_runner.py)."""
    source_content = (REPO_ROOT / env_yaml_relative_path).read_text()
    if SNAKEMAKE_CONDA_DIR.is_dir():
        for marker in SNAKEMAKE_CONDA_DIR.glob("*.yaml"):
            if marker.read_text() == source_content:
                env_hash = marker.stem
                bin_dir = SNAKEMAKE_CONDA_DIR / env_hash / "bin"
                if (SNAKEMAKE_CONDA_DIR / f"{env_hash}.env_setup_done").is_file() and bin_dir.is_dir():
                    return bin_dir
    raise MissingCondaEnvironmentError(f"No se encontro un ambiente Conda ya creado para {env_yaml_relative_path}.")


# Comando usado para consultar la version de cada herramienta externa del
# pipeline, y el archivo de ambiente Conda del que se resuelve (None para
# "snakemake", que se consulta importando el paquete, ver docstring).
#
# "checkm" es otra excepcion: no soporta "--version" (ese flag no es un
# subcomando valido, asi que solo imprime un mensaje de uso/error). La
# version real solo aparece en el encabezado de "checkm -h", por eso su
# comando y su extraccion de version son distintos al resto (ver
# CHECKM_VERSION_PATTERN).
TOOL_VERSION_COMMANDS = {
    "fastp": (["fastp", "--version"], "workflow/envs/fastp.yaml"),
    "spades": (["spades.py", "--version"], "workflow/envs/spades.yaml"),
    "quast": (["quast.py", "--version"], "workflow/envs/quast.yaml"),
    "checkm": (["checkm", "-h"], "workflow/envs/checkm.yaml"),
    "kraken2": (["kraken2", "--version"], "workflow/envs/kraken2.yaml"),
    "prokka": (["prokka", "--version"], "workflow/envs/prokka.yaml"),
    "amrfinder": (["amrfinder", "--version"], "workflow/envs/amrfinder.yaml"),
    "snakemake": (None, None),
}

CHECKM_VERSION_PATTERN = re.compile(r"CheckM v[\d.]+")


def get_tool_version(tool_name: str, command: list[str], env_yaml_relative_path: str) -> str:
    """Ejecuta el comando de version de una herramienta dentro de su propio
    ambiente Conda ya creado. Si ese ambiente todavia no se creo (por
    ejemplo, en un entorno de desarrollo sin todas las herramientas
    instaladas), devuelve "not installed" en vez de fallar: la version es un
    dato de trazabilidad, no debe bloquear el resto del pipeline."""
    try:
        bin_dir = resolve_conda_env_bin(env_yaml_relative_path)
    except MissingCondaEnvironmentError:
        return "not installed"

    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}{os.pathsep}{env.get('PATH', '')}"
    try:
        completed_process = subprocess.run(command, capture_output=True, text=True, check=False, env=env)
    except FileNotFoundError:
        return "not installed"

    # Distintas herramientas imprimen la version en stdout o en stderr.
    version_output = (completed_process.stdout + completed_process.stderr).strip()

    if tool_name == "checkm":
        match = CHECKM_VERSION_PATTERN.search(version_output)
        return match.group(0) if match else "unknown"

    return version_output.splitlines()[0] if version_output else "unknown"


def get_snakemake_version() -> str:
    """Snakemake no vive en un ambiente por-regla: se importa directamente
    (funciona porque este script se invoca con {sys.executable}, el mismo
    interprete que ya tiene Snakemake instalado por definicion, al ser el
    que esta orquestando esta misma corrida)."""
    try:
        import snakemake
        return snakemake.__version__
    except ImportError:
        return "not installed"


def capture_all_tool_versions() -> pd.DataFrame:
    """Consulta la version de cada herramienta del pipeline y las junta en
    una tabla, junto con la fecha en que se consultaron."""
    captured_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    version_rows = []
    for tool_name, (command, env_yaml_relative_path) in TOOL_VERSION_COMMANDS.items():
        if tool_name == "snakemake":
            version = get_snakemake_version()
        else:
            version = get_tool_version(tool_name, command, env_yaml_relative_path)
        version_rows.append({"tool": tool_name, "version": version, "captured_date": captured_date})
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
