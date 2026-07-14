# ============================================================================
# assembly.smk
# Ensamblaje genomico de novo con SPAdes, y filtrado de contigs cortos antes
# de evaluar el ensamblaje (QUAST) o buscar genes de resistencia (AMRFinderPlus).
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
