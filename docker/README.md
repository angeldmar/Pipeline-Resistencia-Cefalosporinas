# Contenedor Docker del pipeline

Empaqueta el pipeline completo —Snakemake, la interfaz web y los ambientes
Conda de todas las herramientas bioinformáticas— en una sola imagen, para
poder usarlo sin instalar nada en el sistema salvo Docker.

## Qué incluye y qué no

| En la imagen | Fuera de la imagen (volumen) |
|---|---|
| Snakemake y el flujo (`workflow/`) | Kraken2 (~8–15 GB) |
| Interfaz web Flask | CheckM (~1.4 GB) |
| Los 12 ambientes Conda ya creados | Datos de entrada (`data/raw/`) |
| Base de datos de AMRFinderPlus (~200 MB) | Resultados (`results/`) |

Las bases grandes (Kraken2, CheckM) y los datos de entrada/salida se montan
como volúmenes: no viajan dentro de la imagen (que sería inmanejable) y
persisten en el host entre corridas. La base de AMRFinderPlus es la
excepción: es pequeña y vive dentro de su propio ambiente Conda (la
herramienta la busca en una ubicación fija, no configurable por ruta), así
que se incluye en la imagen durante el build.

## Requisitos

- Docker con Docker Compose.
- En Mac con chip Apple: la imagen se construye para `linux/amd64` y corre por
  emulación (Docker lo maneja solo). Es más lenta que nativa, pero funciona.
- Espacio en disco: la imagen ocupa **~17 GB** (los 12 ambientes Conda de las
  herramientas) y las bases de datos otros **~17 GB**, es decir del orden de
  **35 GB** en total. Conviene tener ese margen libre antes de empezar.

## Puesta en marcha

Desde la raíz del repositorio:

```bash
# 1. Construir la imagen (crea los ambientes Conda; tarda y pesa varios GB).
docker compose build

# 2. Descargar Kraken2 y CheckM a data/reference (una sola vez, ~17 GB).
#    Cada base se salta si ya esta presente en el volumen.
docker compose run --rm pipeline download-databases

# 3a. Interfaz web: subir un FASTQ/FASTA y ver el reporte.
docker compose up
#     Abrir http://localhost:8080

# 3b. Pipeline por lotes (Snakemake) sobre config/samples.tsv:
docker compose run --rm pipeline pipeline --cores 8
```

## Modos del contenedor

El primer argumento selecciona qué corre (ver `docker/entrypoint.sh`):

| Comando | Qué hace |
|---|---|
| `webapp` (por defecto) | Interfaz web en `0.0.0.0:5000` |
| `pipeline [args]` | `snakemake --use-conda` con los argumentos que sigan |
| `download-databases` | Descarga Kraken2, CheckM y AMRFinderPlus al volumen |
| cualquier otro | Se ejecuta tal cual (`bash`, `pytest`, ...) |

Ejemplos:

```bash
# Un dry-run del pipeline sin ejecutar nada:
docker compose run --rm pipeline pipeline -n --cores 1

# Una shell dentro del contenedor para inspeccionar:
docker compose run --rm pipeline bash

# Correr la suite de pruebas:
docker compose run --rm pipeline pytest -q
```

## Ajustes

- **Índice Kraken2:** por defecto se descarga el `standard-8` (capado a 8 GB).
  Para el índice completo (mejor resolución, más peso) definir `KRAKEN_DB_URL`
  en `docker-compose.yml` o en el entorno antes de `download-databases`.
- **Puerto de la web:** cambiar el mapeo `5000:5000` y la variable
  `WEBAPP_PORT` en `docker-compose.yml`.

## Notas

- Los datos de entrada para la interfaz web se suben por el navegador; para el
  modo por lotes, colocar las lecturas en `data/raw/` y las muestras en
  `config/samples.tsv` (montados desde el host).
- La imagen fija los ambientes Conda al momento del build. Para actualizarlos,
  reconstruir con `docker compose build --no-cache`.
