#!/usr/bin/env bash
#
# Descarga las tres bases de datos de referencia que el pipeline necesita, al
# directorio de referencia (por defecto /pipeline/data/reference, que en uso
# normal es un volumen montado). Se corre una sola vez; cada base se salta si
# ya esta presente, asi que volver a ejecutarlo es seguro y reanudable.
#
#   - CheckM (~1.4 GB): completitud y contaminacion del ensamblaje.
#   - Kraken2 (~8-15 GB segun indice): verificacion taxonomica.
#
# La base de AMRFinderPlus NO se descarga aqui: vive dentro de su ambiente
# Conda (la regla amrfinder usa la ubicacion por defecto de la herramienta,
# no --database) y ya viaja dentro de la imagen desde el build.
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

KRAKEN_DB_URL="${KRAKEN_DB_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz}"
CHECKM_DB_URL="${CHECKM_DB_URL:-https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz}"

mkdir -p "$KRAKEN_DIR" "$CHECKM_DIR"

# ---------------------------------------------------------------------------
# CheckM. Se detecta por un archivo caracteristico de la base (no por un
# marcador propio), para saltar tambien las bases descargadas por otros
# medios que ya esten montadas en el volumen.
# ---------------------------------------------------------------------------
if [ -f "$CHECKM_DIR/taxon_marker_sets.tsv" ]; then
    echo "[CheckM] ya presente, se salta."
else
    echo "[CheckM] descargando (~1.4 GB)..."
    wget -q --show-progress -O /tmp/checkm.tar.gz "$CHECKM_DB_URL"
    tar -xzf /tmp/checkm.tar.gz -C "$CHECKM_DIR"
    rm -f /tmp/checkm.tar.gz
    echo "[CheckM] listo en $CHECKM_DIR"
fi

# ---------------------------------------------------------------------------
# Kraken2 (indice pre-construido; no requiere kraken2-build).
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

echo
echo "Bases de datos (Kraken2, CheckM) listas en $REFERENCE_DIR"
echo "AMRFinderPlus ya viaja dentro de la imagen (ver docker/Dockerfile)."
