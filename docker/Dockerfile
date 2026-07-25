# syntax=docker/dockerfile:1
#
# Imagen del pipeline de deteccion genotipica de resistencia a cefalosporinas
# de tercera generacion en E. coli. Incluye Snakemake, la interfaz web y los
# ambientes Conda de las herramientas bioinformaticas (uno por regla), ya
# creados para linux-64.
#
# Lo que NO va en la imagen (se monta como volumen o se descarga aparte):
#   - Bases de datos de referencia (~17 GB): Kraken2, CheckM, AMRFinderPlus.
#     Se descargan a un volumen con "download-databases" (ver docker/README.md).
#   - Datos de entrada (data/) y resultados (results/): volumenes montados.
#
# Se fija --platform=linux/amd64 porque muchos paquetes de Bioconda solo se
# publican para linux-64/osx-64; en una Mac con chip Apple, Docker lo corre
# por emulacion. Esto sustituye al truco CONDA_SUBDIR=osx-64 que el proyecto
# necesita al instalar en el host (ver README.md).

FROM --platform=linux/amd64 condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="pipeline-amr-ecoli" \
      org.opencontainers.image.description="Deteccion genotipica de resistencia a cefalosporinas de 3G en E. coli" \
      org.opencontainers.image.licenses="MIT"

# Utilidades de sistema que algunas herramientas necesitan en runtime.
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps ca-certificates wget \
    && rm -rf /var/lib/apt/lists/*

# Snakemake, las utilidades base del pipeline y la interfaz web, en el
# ambiente base de Conda (el que queda activo por defecto).
RUN mamba install -y -n base -c conda-forge -c bioconda \
        "snakemake-minimal>=8" pandas pyyaml biopython jinja2 flask \
    && mamba clean -afy

WORKDIR /pipeline

# Copiar el proyecto respetando .dockerignore (sin data/, results/,
# .snakemake/, validacion_estadistica/, etc.).
COPY . /pipeline

# Pre-crear los ambientes Conda de las herramientas (uno por regla), para
# linux-64. Asi la primera corrida no espera a crearlos y la interfaz web los
# encuentra ya listos en .snakemake/conda/. El samples.tsv de ejemplo que
# viene en config/ basta para construir el grafo (DAG) que toca todas las
# reglas; --conda-create-envs-only solo crea los ambientes, no ejecuta nada.
RUN snakemake --use-conda --conda-frontend mamba \
        --conda-create-envs-only --cores 1 \
    && mamba clean -afy

# La anotacion con Prokka no es objetivo por defecto de "all" (ver Snakefile),
# asi que su ambiente no se crea con el paso anterior. Se crea aparte
# apuntando a un output de anotacion, para que la imagen quede completa por si
# alguien corre esa regla. Es un RUN separado para no invalidar la cache de
# los ambientes ya creados arriba al reconstruir.
RUN snakemake --use-conda --conda-frontend mamba \
        --conda-create-envs-only --cores 1 \
        results/tables/annotation/EC001.tsv \
    && mamba clean -afy

# La base de datos de AMRFinderPlus vive DENTRO de su propio ambiente Conda:
# la regla amrfinder no usa --database, sino la ubicacion por defecto de la
# herramienta. Por eso se descarga aqui, en el build, para que viaje en la
# imagen (es pequena, ~200 MB), a diferencia de Kraken2 y CheckM que son
# grandes y van por volumen. Se antepone el bin/ del ambiente al PATH para que
# "amrfinder -u" encuentre blast/hmmer del propio ambiente.
RUN amrfinder_bin_dir="$(dirname "$(find /pipeline/.snakemake/conda -path '*/bin/amrfinder' | head -1)")" \
    && PATH="$amrfinder_bin_dir:$PATH" "$amrfinder_bin_dir/amrfinder" -u

# Puerto de la interfaz web.
EXPOSE 5000

# Entrada/salida y bases de datos van por volumen, nunca dentro de la imagen.
VOLUME ["/pipeline/data", "/pipeline/results"]

ENTRYPOINT ["/pipeline/docker/entrypoint.sh"]
CMD ["webapp"]
