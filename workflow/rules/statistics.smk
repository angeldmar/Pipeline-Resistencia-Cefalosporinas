# ============================================================================
# statistics.smk
# Preparacion de datos (Python) y analisis estadistico formal (R):
# sensibilidad, especificidad, exactitud, kappa, intervalos de confianza y
# coeficiente de variacion. Python solo prepara/limpia; todo el calculo
# estadistico es exclusivo de R (ver run_statistics.R).
# ============================================================================

rule prepare_validation_input:
    # Combina la comparacion contra el estandar de referencia con el
    # desempeno por corrida en el CSV limpio que espera R.
    input:
        reference_comparison="results/tables/reference_comparison.tsv",
        performance_by_sample="results/tables/performance_by_sample.tsv",
    output:
        "results/statistics/validation_input.csv",
    log:
        "logs/prepare_validation_input/prepare.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/prepare_validation_input.py \
          --reference-comparison {input.reference_comparison} \
          --performance-by-sample {input.performance_by_sample} \
          --output {output} \
          > {log} 2>&1
        """


rule run_statistics:
    # Calcula sensibilidad, especificidad, exactitud, kappa, intervalos de
    # confianza y coeficiente de variacion, y genera las graficas finales.
    input:
        "results/statistics/validation_input.csv",
    output:
        confusion_matrix="results/statistics/confusion_matrix.txt",
        classification_metrics="results/statistics/classification_metrics.csv",
        kappa="results/statistics/kappa.csv",
        cv_execution_time="results/statistics/cv_execution_time.csv",
        cv_ram_usage="results/statistics/cv_ram_usage.csv",
    params:
        output_dir="results/statistics",
    log:
        "logs/run_statistics/run.log",
    conda:
        "../envs/r_statistics.yaml"
    shell:
        """
        Rscript workflow/scripts/run_statistics.R {input} {params.output_dir} \
          > {log} 2>&1
        """
