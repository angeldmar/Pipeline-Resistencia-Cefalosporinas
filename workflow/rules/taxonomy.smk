# ============================================================================
# taxonomy.smk
# Identificacion taxonomica con Kraken2, para confirmar que cada muestra
# corresponde a Escherichia coli antes de reportar genes de resistencia.
# ============================================================================

rule kraken2:
    # Clasifica taxonomicamente las lecturas ya recortadas/filtradas contra
    # una base de datos de referencia. Se corre sobre lecturas (no sobre el
    # ensamblaje) porque es el uso estandar de Kraken2 y detecta mezclas de
    # especies que un ensamblaje ya colapsado podria enmascarar.
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
    output:
        report="results/taxonomy/kraken2/{sample}/report.tsv",
        classification="results/taxonomy/kraken2/{sample}/classification.tsv",
    params:
        database=config["paths"]["kraken_database"],
    log:
        "logs/kraken2/{sample}.log",
    conda:
        "../envs/kraken2.yaml"
    threads:
        config["threads"]["kraken2"]
    shell:
        """
        kraken2 \
          --db {params.database} \
          --paired {input.r1} {input.r2} \
          --threads {threads} \
          --report {output.report} \
          --output {output.classification} \
          > {log} 2>&1
        """


rule parse_kraken2:
    # Extrae del reporte de Kraken2 el % asignado a E. coli, el taxon
    # predominante y el % de otras especies, y clasifica la muestra en
    # PASS/WARNING/FAIL. Las lecturas asignadas a Shigella se excluyen del
    # calculo de contaminacion (ver nota en parse_kraken2.py) y en su lugar
    # marcan la muestra para revision manual, por la cercania genomica entre
    # ambos generos.
    #
    # Igual que en los parsers anteriores, se escribe un archivo POR MUESTRA
    # para evitar condiciones de carrera cuando varias muestras corren en
    # paralelo bajo Snakemake.
    input:
        "results/taxonomy/kraken2/{sample}/report.tsv",
    output:
        "results/tables/taxonomy/{sample}.tsv",
    log:
        "logs/parse_kraken2/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_kraken2.py parse {wildcards.sample} {input} \
          --output-dir results/tables/taxonomy \
          --minimum-ecoli-percentage {config[taxonomy][minimum_ecoli_percentage]} \
          --warning-ecoli-percentage {config[taxonomy][warning_ecoli_percentage]} \
          --maximum-contaminant-percentage {config[taxonomy][maximum_contaminant_percentage]} \
          > {log} 2>&1
        """
