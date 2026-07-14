# ============================================================================
# quality_control.smk
# Control de calidad de lecturas crudas con fastp: recorte de adaptadores,
# filtrado por calidad/longitud y generacion de reportes HTML/JSON.
# ============================================================================

rule fastp:
    # Ejecuta fastp sobre las lecturas crudas paired-end de una muestra.
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        json="results/qc/fastp/{sample}.json",
        html="results/qc/fastp/{sample}.html",
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yaml"
    threads:
        config["threads"]["fastp"]
    shell:
        """
        fastp \
          --in1 {input.r1} \
          --in2 {input.r2} \
          --out1 {output.r1} \
          --out2 {output.r2} \
          --json {output.json} \
          --html {output.html} \
          --length_required {config[quality][minimum_length]} \
          --qualified_quality_phred {config[quality][minimum_phred]} \
          --thread {threads} \
          > {log} 2>&1
        """


rule parse_fastp:
    # Extrae del JSON de fastp las metricas que le interesan al pipeline
    # (lecturas retenidas, bases Q20/Q30, GC, duplicacion, % filtrado) y
    # estima la cobertura del genoma a partir de esas mismas lecturas,
    # clasificandola en PASS/WARNING/FAIL segun los umbrales de config.yaml.
    # Todo queda en una tabla individual por muestra.
    #
    # Se escribe un archivo POR MUESTRA (no se acumula directamente en un
    # unico fastp_summary.tsv) porque varias muestras pueden ejecutarse en
    # paralelo bajo Snakemake, y varios procesos escribiendo a la vez sobre
    # el mismo archivo compartido causaria condiciones de carrera. La tabla
    # combinada de todas las muestras se arma en un paso de agregacion aparte
    # (ver seccion de integracion de resultados).
    input:
        json="results/qc/fastp/{sample}.json",
    output:
        "results/tables/fastp/{sample}.tsv",
    log:
        "logs/parse_fastp/{sample}.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/parse_fastp.py parse {wildcards.sample} {input.json} \
          --output-dir results/tables/fastp \
          --genome-size {config[quality][estimated_genome_size]} \
          --minimum-coverage {config[quality][minimum_coverage]} \
          --warning-coverage {config[quality][warning_coverage]} \
          > {log} 2>&1
        """
