#!/usr/bin/env bash
#
# Descarga las tres bases de datos de referencia que el pipeline necesita, al
# directorio de referencia (por defecto /pipeline/data/reference, que en uso
# normal es un volumen montado). Se corre una sola vez; cada base se salta si
# ya esta presente, asi que volver a ejecutarlo es seguro y reanudable.
#
#   - CheckM (~1.4 GB): completitud y contaminacion del ensamblaje.
#   - Kraken2 (~8-15 GB segun indice): verificacion taxonomica.
#   - AMRFinderPlus (~200 MB): deteccion de genes de resistencia.
#
# Variables de entorno para ajustar:
#   REFERENCE_DIR   destino (por defecto /pipeline/data/reference)
#   KRAKEN_DB_URL   indice Kraken2 a descargar. Por defecto el "standard-8"
#                   (capado a 8 GB); para el completo (~mas resolucion, mas
#                   peso) usar el indice k2_standard sin el sufijo _08gb.

set -euo pipefail

REFERENCE_DIR="${REFERENCE_DIR:-/pipeline/data/reference}"
KRAKEN_DIR="$REFERENCE_DIR/kraken2"
CHECKM_DIR="$REFERENCE_DIR/checkm"
AMRFINDER_DIR="$REFERENCE_DIR/amrfinder"

KRAKEN_DB_URL="${KRAKEN_DB_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz}"
CHECKM_DB_URL="${CHECKM_DB_URL:-https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz}"

mkdir -p "$KRAKEN_DIR" "$CHECKM_DIR" "$AMRFINDER_DIR"

# ---------------------------------------------------------------------------
# CheckM
# ---------------------------------------------------------------------------
if [ -f "$CHECKM_DIR/.download_done" ]; then
    echo "[CheckM] ya presente, se salta."
else
    echo "[CheckM] descargando (~1.4 GB)..."
    wget -q --show-progress -O /tmp/checkm.tar.gz "$CHECKM_DB_URL"
    tar -xzf /tmp/checkm.tar.gz -C "$CHECKM_DIR"
    rm -f /tmp/checkm.tar.gz
    touch "$CHECKM_DIR/.download_done"
    echo "[CheckM] listo en $CHECKM_DIR"
fi

# ---------------------------------------------------------------------------
# Kraken2 (indice pre-construido; no requiere kraken2-build)
# ---------------------------------------------------------------------------
if [ -f "$KRAKEN_DIR/hash.k2d" ]; then
    echo "[Kraken2] ya presente, se salta."
else
    echo "[Kraken2] descargando indice desde $KRAKEN_DB_URL ..."
    wget -q --show-progress -O /tmp/kraken2.tar.gz "$KRAKEN_DB_URL"
    tar -xzf /tmp/kraken2.tar.gz -C "$KRAKEN_DIR"
    rm -f /tmp/kraken2.tar.gz
    echo "[Kraken2] listo en $KRAKEN_DIR"
fi

# ---------------------------------------------------------------------------
# AMRFinderPlus (se descarga con la herramienta del ambiente Conda ya creado)
# ---------------------------------------------------------------------------
if [ -d "$AMRFINDER_DIR/latest" ]; then
    echo "[AMRFinderPlus] ya presente, se salta."
else
    echo "[AMRFinderPlus] descargando base de datos..."
    amrfinder_update_bin=$(find /pipeline/.snakemake/conda -path "*/bin/amrfinder_update" 2>/dev/null | head -1)
    if [ -z "$amrfinder_update_bin" ]; then
        echo "ERROR: no se encontro amrfinder_update en los ambientes Conda." >&2
        echo "       Construya la imagen con los ambientes pre-creados (ver docker/README.md)." >&2
        exit 1
    fi
    "$amrfinder_update_bin" --database "$AMRFINDER_DIR"
    echo "[AMRFinderPlus] listo en $AMRFINDER_DIR"
fi

echo
echo "Bases de datos listas en $REFERENCE_DIR"
