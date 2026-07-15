# ============================================================================
# amr_detection.smk
# Deteccion de genes de resistencia antimicrobiana con AMRFinderPlus, y
# clasificacion mecanistica de las beta-lactamasas detectadas (BLEE, AmpC,
# carbapenemasa) segun config/resistance_targets.yaml.
# ============================================================================

rule amrfinder:
    # Busca genes de resistencia directamente sobre el ensamblaje de
    # nucleotidos filtrado; no depende de la anotacion de Prokka.
    input:
        assembly="results/assemblies/{sample}/contigs.filtered.fasta",
    output:
        table="results/amr/amrfinder/{sample}.tsv",
    log:
        "logs/amrfinder/{sample}.log",
    conda:
        "../envs/amrfinder.yaml"
    threads:
        config["threads"]["amrfinder"]
    shell:
        """
        amrfinder \
          --nucleotide {input.assembly} \
          --organism Escherichia \
          --threads {threads} \
          --output {output.table} \
          > {log} 2>&1
        """


rule parse_amrfinder:
    # Normaliza la salida cruda de AMRFinderPlus a un esquema de columnas
    # propio y estable, marcando (sin descartar) las detecciones que no
    # alcanzan los umbrales minimos de identidad/cobertura de config.yaml.
    #
    # Igual que en los demas parsers, se escribe un archivo POR MUESTRA para
    # evitar condiciones de carrera cuando varias muestras corren en paralelo.
    input:
        "results/amr/amrfinder/{sample}.tsv",
    output:
        "results/tables/amr/{sample}.tsv",
    log:
        "logs/parse_amrfinder/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_amrfinder.py parse {wildcards.sample} {input} \
          --output-dir results/tables/amr \
          --minimum-identity {config[amr][minimum_identity]} \
          --minimum-gene-coverage {config[amr][minimum_gene_coverage]} \
          > {log} 2>&1
        """


rule classify_cephalosporin_genes:
    # Paso de agregacion final (no por muestra): toma el listado largo ya
    # combinado de TODAS las muestras y le agrega la categoria mecanistica de
    # cada beta-lactamasa detectada (BLEE/AmpC/carbapenemasa/Other).
    input:
        amr_summary="results/tables/amr_summary.tsv",
        resistance_targets=config["resistance_targets"],
    output:
        "results/tables/amr_classified.tsv",
    log:
        "logs/classify_cephalosporin_genes/classify.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/classify_cephalosporin_genes.py {input.amr_summary} \
          --resistance-targets {input.resistance_targets} \
          --output {output} \
          > {log} 2>&1
        """
