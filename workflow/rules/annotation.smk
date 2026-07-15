# ============================================================================
# annotation.smk
# Anotacion genomica con Prokka. Es un paso informativo/complementario: la
# deteccion de AMR (amr_detection.smk) usa AMRFinderPlus directamente sobre
# el ensamblaje de nucleotidos y NO depende de esta anotacion.
# ============================================================================

rule prokka:
    # Anota el ensamblaje filtrado: genes codificantes (CDS), ARN ribosomal
    # y de transferencia, etc. --force permite re-ejecutar sobre una carpeta
    # de salida existente (por ejemplo, al reprocesar una muestra).
    input:
        "results/assemblies/{sample}/contigs.filtered.fasta",
    output:
        gff="results/annotation/{sample}/{sample}.gff",
        gbk="results/annotation/{sample}/{sample}.gbk",
        faa="results/annotation/{sample}/{sample}.faa",
        ffn="results/annotation/{sample}/{sample}.ffn",
        tsv="results/annotation/{sample}/{sample}.tsv",
        txt="results/annotation/{sample}/{sample}.txt",
    params:
        outdir="results/annotation/{sample}",
        prefix="{sample}",
    log:
        "logs/prokka/{sample}.log",
    conda:
        "../envs/prokka.yaml"
    threads:
        config["threads"]["prokka"]
    shell:
        """
        prokka \
          --outdir {params.outdir} \
          --prefix {params.prefix} \
          --genus Escherichia \
          --species coli \
          --strain {wildcards.sample} \
          --cpus {threads} \
          --force \
          {input} \
          > {log} 2>&1
        """


rule parse_prokka:
    # Organiza los resultados de Prokka en una tabla de metricas por muestra:
    # numero de CDS, ARN ribosomal, ARN de transferencia, genes hipoteticos
    # (sin funcion conocida) y version de la herramienta usada.
    #
    # Igual que en los demas parsers, se escribe un archivo POR MUESTRA para
    # evitar condiciones de carrera cuando varias muestras corren en paralelo.
    input:
        summary="results/annotation/{sample}/{sample}.txt",
        annotation_tsv="results/annotation/{sample}/{sample}.tsv",
    output:
        "results/tables/annotation/{sample}.tsv",
    log:
        "logs/parse_prokka/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_prokka.py parse {wildcards.sample} \
          --summary-txt {input.summary} \
          --annotation-tsv {input.annotation_tsv} \
          --output-dir results/tables/annotation \
          > {log} 2>&1
        """
