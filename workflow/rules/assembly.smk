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
    # --phred-offset 33 fijo: sin este flag, SPAdes intenta autodetectar el
    # offset a partir de la variabilidad de calidades en las lecturas, y
    # falla ("Failed to determine offset!") si esa variabilidad es baja.
    # Phred+33 es el estandar universal en datos Illumina desde ~2011 (el
    # antiguo Phred+64 esta obsoleto), asi que fijarlo evita depender de esa
    # heuristica sin perder generalidad.
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
    output:
        assembly="results/assemblies/{sample}/contigs.fasta",
        performance="results/tables/performance/{sample}_spades.tsv",
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
          --output {output.performance} \
          -- \
          spades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {params.outdir} \
          -t {threads} \
          --careful \
          --phred-offset 33 \
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
        performance="results/tables/performance/{sample}_quast.tsv",
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
          --output {output.performance} \
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


rule combine_quast:
    # Junta las tablas individuales de TODAS las muestras en un unico resumen.
    input:
        expand("results/tables/quast/{sample}.tsv", sample=SAMPLES),
    output:
        "results/tables/quast_summary.tsv",
    log:
        "logs/combine_quast/combine.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_quast.py combine \
          --input-dir results/tables/quast \
          --output {output} \
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
        performance="results/tables/performance/{sample}_checkm.tsv",
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
        # {sys.executable} en vez de "python": el ambiente de checkm no
        # puede tener un Python moderno (ver comentario en
        # workflow/envs/checkm.yaml), asi que run_with_timing.py se invoca
        # con el interprete que lanzo Snakemake -- que no depende del
        # ambiente activado para esta regla -- en vez de con el Python 2.7
        # que trae CheckM como dependencia interna.
        """
        mkdir -p {params.bin_dir}
        cp {input} {params.bin_dir}/{wildcards.sample}.fasta
        checkm data setRoot {params.database} > {log} 2>&1
        {sys.executable} workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module checkm \
          --threads {threads} \
          --output {output.performance} \
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


rule combine_checkm:
    # Junta las tablas individuales de TODAS las muestras en un resumen, y
    # ademas genera el registro de exclusiones (muestras FAIL, ver seccion
    # 10 del diseno del pipeline: nunca se descartan en silencio).
    input:
        expand("results/tables/checkm/{sample}.tsv", sample=SAMPLES),
    output:
        summary="results/tables/checkm_summary.tsv",
        exclusions="results/tables/checkm_exclusions.tsv",
    log:
        "logs/combine_checkm/combine.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_checkm.py combine \
          --input-dir results/tables/checkm \
          --output {output.summary} \
          --exclusions-output {output.exclusions} \
          > {log} 2>&1
        """
