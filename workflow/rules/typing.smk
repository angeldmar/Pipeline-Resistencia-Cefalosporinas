# ============================================================================
# typing.smk
# Tipificacion de secuencia multilocus (MLST): da contexto epidemiologico
# (linaje/secuencia tipo, ej. ST131) a los genes de resistencia detectados.
# No interactua con la deteccion de AMR; es un modulo independiente que
# enriquece la tabla maestra con una columna de contexto adicional.
# ============================================================================

rule mlst:
    # mlst imprime su resultado por stdout (no tiene flag --output), por eso
    # la regla redirige stdout directamente al archivo de resultados (igual
    # que la regla abricate).
    input:
        assembly="results/assemblies/{sample}/contigs.filtered.fasta",
    output:
        table="results/typing/mlst/{sample}.tsv",
        performance="results/tables/performance/{sample}_mlst.tsv",
    log:
        "logs/mlst/{sample}.log",
    conda:
        "../envs/mlst.yaml"
    threads:
        config["threads"]["mlst"]
    shell:
        """
        python workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module mlst \
          --threads {threads} \
          --output {output.performance} \
          -- \
          mlst \
          --scheme {config[mlst][scheme]} \
          --threads {threads} \
          {input.assembly} \
          > {output.table} 2> {log}
        """


rule parse_mlst:
    # Normaliza la fila cruda de mlst (sin encabezado, columnas de alelo
    # variables segun el esquema) al esquema del pipeline.
    input:
        "results/typing/mlst/{sample}.tsv",
    output:
        "results/tables/mlst/{sample}.tsv",
    log:
        "logs/parse_mlst/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_mlst.py parse {wildcards.sample} {input} \
          --output-dir results/tables/mlst \
          > {log} 2>&1
        """


rule combine_mlst:
    # Junta las tablas individuales de TODAS las muestras en un unico resumen.
    input:
        expand("results/tables/mlst/{sample}.tsv", sample=SAMPLES),
    output:
        "results/tables/mlst_summary.tsv",
    log:
        "logs/combine_mlst/combine.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_mlst.py combine \
          --input-dir results/tables/mlst \
          --output {output} \
          > {log} 2>&1
        """
