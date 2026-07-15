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

rule merge_results:
    # Combina los resultados de todos los modulos en una tabla maestra (una
    # fila por muestra). No mete el detalle gen-por-gen en columnas: eso
    # sigue viviendo en la tabla larga de AMR (amr_summary.tsv).
    input:
        samples=config["samples"],
        fastp_summary="results/tables/fastp_summary.tsv",
        quast_summary="results/tables/quast_summary.tsv",
        checkm_summary="results/tables/checkm_summary.tsv",
        taxonomy_summary="results/tables/taxonomy_summary.tsv",
        amr_summary="results/tables/amr_summary.tsv",
        reference_comparison="results/tables/reference_comparison.tsv",
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
          --output {output} \
          > {log} 2>&1
        """
