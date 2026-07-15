# ============================================================================
# reports.smk
# Integracion final de resultados (tabla maestra) y generacion de reportes.
#
# NOTA: al igual que las reglas de agregacion en amr_detection.smk, esta
# regla depende de los *_summary.tsv de cada modulo (fastp, QUAST, CheckM,
# Kraken2, AMRFinderPlus), que a su vez se arman juntando la tabla de TODAS
# las muestras. Esa agregacion por muestra completa necesita SAMPLES, que
# recien queda definido al construir el Snakefile principal. Por ahora,
# master_results.tsv se genera manualmente con merge_results.py y esta regla
# se conectara al grafo completo cuando se arme el Snakefile principal.
# ============================================================================

rule combine_performance:
    # Junta los registros de tiempo/CPU/RAM de cada ejecucion individual
    # (results/tables/performance/{sample}_{module}.tsv, generados por
    # run_with_timing.py) en un listado largo y dos resumenes derivados
    # (por muestra y por modulo). Ver combine_performance.py para el detalle.
    input:
        performance_dir="results/tables/performance",
    output:
        summary="results/tables/performance_summary.tsv",
        by_sample="results/tables/performance_by_sample.tsv",
        by_module="results/tables/performance_by_module.tsv",
    log:
        "logs/combine_performance/combine.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/combine_performance.py \
          --input-dir {input.performance_dir} \
          --output {output.summary} \
          --by-sample-output {output.by_sample} \
          --by-module-output {output.by_module} \
          > {log} 2>&1
        """


rule merge_results:
    # Combina los resultados de todos los modulos en una tabla maestra (una
    # fila por muestra). No mete el detalle gen-por-gen en columnas: eso
    # sigue viviendo en la tabla larga de AMR (amr_summary.tsv). Tampoco mete
    # el desglose de tiempo/RAM por modulo: solo el total por muestra (ver
    # combine_performance.py).
    input:
        samples=config["samples"],
        fastp_summary="results/tables/fastp_summary.tsv",
        quast_summary="results/tables/quast_summary.tsv",
        checkm_summary="results/tables/checkm_summary.tsv",
        taxonomy_summary="results/tables/taxonomy_summary.tsv",
        amr_summary="results/tables/amr_summary.tsv",
        reference_comparison="results/tables/reference_comparison.tsv",
        performance_by_sample="results/tables/performance_by_sample.tsv",
    output:
        "results/tables/master_results.tsv",
    log:
        "logs/merge_results/merge.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/merge_results.py \
          --samples {input.samples} \
          --fastp-summary {input.fastp_summary} \
          --quast-summary {input.quast_summary} \
          --checkm-summary {input.checkm_summary} \
          --taxonomy-summary {input.taxonomy_summary} \
          --amr-summary {input.amr_summary} \
          --reference-comparison {input.reference_comparison} \
          --performance-by-sample {input.performance_by_sample} \
          --output {output} \
          > {log} 2>&1
        """


rule capture_tool_versions:
    # Registra, UNA sola vez por corrida (no por muestra: la version de una
    # herramienta no cambia entre muestras), las versiones instaladas de
    # cada herramienta externa del pipeline.
    output:
        "data/metadata/tool_versions.tsv",
    log:
        "logs/capture_tool_versions/capture.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/capture_tool_versions.py --output {output} > {log} 2>&1
        """


rule generate_report:
    # Genera el reporte HTML individual de una muestra: identificador,
    # calidad, cobertura, ensamblaje, taxonomia, completitud/contaminacion,
    # genes de resistencia con su interpretacion GENOTIPICA (nunca una
    # conclusion clinica), comparacion con la referencia, desempeno y
    # versiones de herramientas. Ver generate_report.py.
    input:
        master_results="results/tables/master_results.tsv",
        amr_classified="results/tables/amr_classified.tsv",
        performance_summary="results/tables/performance_summary.tsv",
        tool_versions="data/metadata/tool_versions.tsv",
    output:
        "results/reports/{sample}.html",
    log:
        "logs/generate_report/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/generate_report.py \
          --sample-id {wildcards.sample} \
          --master-results {input.master_results} \
          --amr-classified {input.amr_classified} \
          --performance-summary {input.performance_summary} \
          --tool-versions {input.tool_versions} \
          --output {output} \
          > {log} 2>&1
        """
