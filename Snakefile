# ============================================================================
# Snakefile
# Punto de entrada del pipeline de deteccion de resistencia a cefalosporinas
# de tercera generacion en E. coli. Orquesta, en orden, todos los modulos ya
# construidos y probados individualmente (ver README.md para el detalle de
# cada parte): validacion de metadatos, descarga, control de calidad,
# ensamblaje, taxonomia, anotacion, deteccion de AMR, integracion de
# resultados, estadistica en R y generacion de reportes.
# ============================================================================

import sys
from pathlib import Path

import pandas as pd

configfile: "config/config.yaml"

# ----------------------------------------------------------------------------
# Validacion de metadatos y lista de muestras
# ----------------------------------------------------------------------------
# Se reutiliza directamente la funcion validate_samples() (la misma que usa
# el script de linea de comandos de la parte 3), en vez de duplicar su
# logica aqui. Si samples.tsv tiene algun problema (columnas faltantes,
# duplicados, fenotipos invalidos, accesiones mal formadas, etc.), Snakemake
# falla inmediatamente al cargar el Snakefile, antes de intentar construir
# ningun grafo de dependencias con datos no confiables.
sys.path.insert(0, str(Path("workflow/scripts")))
from validate_samples import validate_samples

samples_df = validate_samples(Path(config["samples"]))
SAMPLES = samples_df["sample_id"].tolist()

# ----------------------------------------------------------------------------
# Reglas de cada modulo. SAMPLES debe existir ANTES de este punto: varias
# reglas de agregacion (combine_fastp, combine_quast, combine_checkm,
# combine_taxonomy, combine_amr, combine_abricate) usan expand(..., sample=SAMPLES)
# para depender de la tabla de TODAS las muestras a la vez.
# ----------------------------------------------------------------------------
include: "workflow/rules/download.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/amr_detection.smk"
include: "workflow/rules/typing.smk"
include: "workflow/rules/statistics.smk"
include: "workflow/rules/reports.smk"

# ----------------------------------------------------------------------------
# Objetivo final: un reporte HTML por muestra, la tabla maestra, los
# registros de exclusion/revision manual (para que ninguna muestra excluida
# desaparezca en silencio), y los resultados del analisis estadistico en R.
#
# La anotacion (Prokka) queda disponible como regla (resultados/tables/
# annotation/{sample}.tsv) pero NO es un objetivo por defecto de "all": es
# informativa/complementaria (AMRFinderPlus no depende de ella) y no forma
# parte de los campos que exige el reporte individual (seccion 19).
# ----------------------------------------------------------------------------
rule all:
    input:
        expand("results/reports/{sample}.html", sample=SAMPLES),
        "results/tables/master_results.tsv",
        "results/tables/amr_classified.tsv",
        "results/tables/checkm_exclusions.tsv",
        "results/tables/taxonomy_manual_review.tsv",
        "results/tables/engine_concordance.tsv",
        "results/tables/mlst_summary.tsv",
        "results/statistics/confusion_matrix.txt",
        "results/statistics/classification_metrics.csv",
        "results/statistics/kappa.csv",
        "results/statistics/cv_execution_time.csv",
        "results/statistics/cv_ram_usage.csv",
        "results/statistics/engine_concordance_kappa.csv",
