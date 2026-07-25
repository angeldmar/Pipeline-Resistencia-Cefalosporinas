#!/usr/bin/env bash
#
# Punto de entrada del contenedor. Selecciona que se ejecuta segun el primer
# argumento; el ambiente base de Conda (con Snakemake y la interfaz web) ya
# esta activo en la imagen.
#
#   webapp             -> interfaz web Flask en 0.0.0.0:5000 (por defecto)
#   pipeline [args]    -> Snakemake con --use-conda y los argumentos que sigan
#   download-databases -> descarga Kraken2, CheckM y AMRFinderPlus al volumen
#   <cualquier otro>   -> se ejecuta tal cual (bash, snakemake, pytest, ...)

set -euo pipefail

# El env base ya trae Snakemake, Flask y utilidades; activarlo hace que
# "conda"/"snakemake" esten en el PATH para todos los modos.
source /opt/conda/etc/profile.d/conda.sh
conda activate base

command="${1:-webapp}"

case "$command" in
    webapp)
        export WEBAPP_HOST="${WEBAPP_HOST:-0.0.0.0}"
        export WEBAPP_PORT="${WEBAPP_PORT:-5000}"
        export WEBAPP_DEBUG="${WEBAPP_DEBUG:-0}"
        echo "Interfaz web disponible en http://localhost:${WEBAPP_PORT}"
        exec python webapp/app.py
        ;;
    pipeline)
        shift
        # Valores por defecto razonables; se pueden sobreescribir pasando mas
        # argumentos (p. ej. "pipeline --cores 8 --config samples=...").
        exec snakemake --use-conda --conda-frontend mamba "$@"
        ;;
    download-databases)
        shift
        exec /pipeline/docker/download_databases.sh "$@"
        ;;
    *)
        # Cualquier otro comando se ejecuta tal cual dentro del env base.
        exec "$@"
        ;;
esac
