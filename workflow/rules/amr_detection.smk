# ============================================================================
# amr_detection.smk
# Deteccion de genes de resistencia antimicrobiana con AMRFinderPlus,
# clasificacion mecanistica de las beta-lactamasas detectadas (BLEE, AmpC,
# carbapenemasa) segun config/resistance_targets.yaml, y comparacion contra
# el estandar de referencia documentado en samples.tsv.
#
# NOTA: las reglas de agregacion (classify_cephalosporin_genes y
# compare_to_reference) toman como entrada results/tables/amr_summary.tsv,
# que se arma juntando la tabla de CADA muestra (parse_amrfinder.py combine).
# Esa regla de "combinar todas las muestras" necesita conocer la lista
# completa de muestras (SAMPLES), que recien queda definida al construir el
# Snakefile principal (parte pendiente del roadmap). Por eso, por ahora,
# amr_summary.tsv se genera manualmente con:
#   python workflow/scripts/parse_amrfinder.py combine
# y se conectara al grafo de Snakemake cuando se arme el Snakefile principal.
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


rule compare_to_reference:
    # Paso de agregacion final: compara, muestra por muestra, el gen
    # esperado (samples.tsv, columna expected_genes) contra los genes que
    # realmente detecto el pipeline, y arma la matriz TP/TN/FP/FN que R usara
    # para calcular sensibilidad, especificidad y kappa. Ningun calculo
    # estadistico se hace aqui, solo se prepara y clasifica la matriz de datos.
    input:
        samples=config["samples"],
        amr_summary="results/tables/amr_summary.tsv",
    output:
        "results/tables/reference_comparison.tsv",
    log:
        "logs/compare_to_reference/compare.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/compare_to_reference.py \
          --samples {input.samples} \
          --amr-table {input.amr_summary} \
          --output {output} \
          > {log} 2>&1
        """
