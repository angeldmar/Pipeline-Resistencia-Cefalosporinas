# ============================================================================
# assembly.smk
# Ensamblaje genomico de novo con SPAdes, filtrado de contigs cortos, y
# evaluacion del ensamblaje con QUAST (numero de contigs, N50, longitud
# total, contenido GC) con clasificacion automatica PASS/WARNING/FAIL.
# ============================================================================

rule spades:
    # Ensambla de novo las lecturas ya recortadas/filtradas (results/trimmed/)
    # en contigs. --careful reduce errores de ensamblaje (mismatches/indels)
    # a costa de mayor tiempo de computo, adecuado para genomas bacterianos.
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
    output:
        assembly="results/assemblies/{sample}/contigs.fasta",
    params:
        outdir="results/assemblies/{sample}",
    log:
        "logs/spades/{sample}.log",
    conda:
        "../envs/spades.yaml"
    threads:
        config["threads"]["spades"]
    shell:
        """
        spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {params.outdir} \
          -t {threads} \
          --careful \
          > {log} 2>&1
        """


rule filter_contigs:
    # Descarta los contigs mas cortos que assembly.minimum_contig_length
    # (config.yaml), ya que suelen ser ruido de ensamblaje mas que secuencia
    # genomica confiable. La salida (contigs.filtered.fasta) es la que usan
    # QUAST y AMRFinderPlus mas adelante.
    input:
        "results/assemblies/{sample}/contigs.fasta",
    output:
        "results/assemblies/{sample}/contigs.filtered.fasta",
    log:
        "logs/filter_contigs/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/filter_contigs.py {input} {output} \
          --minimum-length {config[assembly][minimum_contig_length]} \
          > {log} 2>&1
        """


rule quast:
    # Evalua el ensamblaje filtrado con QUAST: numero de contigs, N50,
    # longitud total, contig mas largo y contenido GC.
    input:
        "results/assemblies/{sample}/contigs.filtered.fasta",
    output:
        report="results/qc/quast/{sample}/report.tsv",
    params:
        outdir="results/qc/quast/{sample}",
    log:
        "logs/quast/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads:
        config["threads"]["quast"]
    shell:
        """
        quast.py {input} \
          --output-dir {params.outdir} \
          --threads {threads} \
          > {log} 2>&1
        """


rule parse_quast:
    # Extrae del report.tsv de QUAST las metricas clave del ensamblaje y lo
    # clasifica en PASS/WARNING/FAIL segun los umbrales de config.yaml.
    #
    # Igual que en parse_fastp, se escribe un archivo POR MUESTRA (no se
    # acumula directamente en un unico quast_summary.tsv) para evitar
    # condiciones de carrera cuando varias muestras corren en paralelo. La
    # tabla combinada se arma en un paso de agregacion aparte.
    input:
        "results/qc/quast/{sample}/report.tsv",
    output:
        "results/tables/quast/{sample}.tsv",
    log:
        "logs/parse_quast/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_quast.py parse {wildcards.sample} {input} \
          --output-dir results/tables/quast \
          --maximum-contigs {config[assembly][maximum_contigs]} \
          --minimum-total-length {config[assembly][minimum_total_length]} \
          --maximum-total-length {config[assembly][maximum_total_length]} \
          --n50-warning-threshold {config[assembly][n50_warning_threshold]} \
          > {log} 2>&1
        """
