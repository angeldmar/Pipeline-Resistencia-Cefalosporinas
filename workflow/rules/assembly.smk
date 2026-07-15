# ============================================================================
# assembly.smk
# Ensamblaje genomico de novo con SPAdes, filtrado de contigs cortos,
# evaluacion del ensamblaje con QUAST (numero de contigs, N50, longitud
# total, contenido GC), y estimacion de completitud/contaminacion con CheckM.
# Todas estas reglas evaluan la calidad del mismo ensamblaje filtrado, por
# eso se mantienen juntas en este archivo.
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
        python workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module spades \
          --threads {threads} \
          --output results/tables/performance/{wildcards.sample}_spades.tsv \
          -- \
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
        python workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module quast \
          --threads {threads} \
          --output results/tables/performance/{wildcards.sample}_quast.tsv \
          -- \
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


rule checkm:
    # Estima completitud y contaminacion con CheckM, a partir de genes
    # marcadores de copia unica especificos del linaje taxonomico. CheckM
    # espera una CARPETA de bins (genomas), no un solo archivo, asi que la
    # regla primero copia el ensamblaje filtrado de esta muestra a su propia
    # carpeta de bins antes de correr "checkm lineage_wf".
    input:
        "results/assemblies/{sample}/contigs.filtered.fasta",
    output:
        report="results/qc/checkm/{sample}/checkm_summary.tsv",
    params:
        bin_dir="results/qc/checkm/{sample}/bins",
        outdir="results/qc/checkm/{sample}/output",
        database=config["paths"]["checkm_database"],
    log:
        "logs/checkm/{sample}.log",
    conda:
        "../envs/checkm.yaml"
    threads:
        config["threads"]["checkm"]
    shell:
        """
        mkdir -p {params.bin_dir}
        cp {input} {params.bin_dir}/{wildcards.sample}.fasta
        checkm data setRoot {params.database} > {log} 2>&1
        python workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module checkm \
          --threads {threads} \
          --output results/tables/performance/{wildcards.sample}_checkm.tsv \
          -- \
          checkm lineage_wf \
          -x fasta \
          --tab_table \
          -f {output.report} \
          -t {threads} \
          {params.bin_dir} {params.outdir} \
          >> {log} 2>&1
        """


rule parse_checkm:
    # Extrae completitud y contaminacion del reporte de CheckM y clasifica la
    # muestra en PASS/FAIL segun los umbrales de config.yaml. Sigue el mismo
    # patron de archivo-por-muestra que parse_fastp y parse_quast.
    input:
        "results/qc/checkm/{sample}/checkm_summary.tsv",
    output:
        "results/tables/checkm/{sample}.tsv",
    log:
        "logs/parse_checkm/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_checkm.py parse {wildcards.sample} {input} \
          --output-dir results/tables/checkm \
          --minimum-completeness {config[assembly][minimum_completeness]} \
          --maximum-contamination {config[assembly][maximum_contamination]} \
          > {log} 2>&1
        """
