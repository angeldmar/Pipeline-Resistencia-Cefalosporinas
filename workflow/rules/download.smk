# ============================================================================
# download.smk
# Descarga de lecturas crudas desde SRA (ver workflow/scripts/download_data.py).
# Envuelta con run_with_timing.py para registrar tiempo/CPU/RAM de la descarga,
# igual que el resto de los modulos "importantes" del pipeline.
# ============================================================================

rule download_sample:
    # Descarga, comprime y verifica las lecturas paired-end de UNA muestra.
    # download_data.py ya soporta --sample-id para procesar una sola muestra
    # (ver parte de descarga de datos), lo que encaja naturalmente con el
    # modelo de reglas-por-muestra de Snakemake.
    input:
        samples=config["samples"],
    output:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz",
        performance="results/tables/performance/{sample}_download.tsv",
    log:
        "logs/download/{sample}.log",
    conda:
        "../envs/python.yaml"
    threads:
        config["threads"]["download"]
    shell:
        """
        python workflow/scripts/run_with_timing.py \
          --sample-id {wildcards.sample} \
          --module download \
          --threads {threads} \
          --output {output.performance} \
          -- \
          python workflow/scripts/download_data.py {input.samples} \
            --sample-id {wildcards.sample} \
            --raw-dir data/raw \
            --manifest data/metadata/download_manifest.tsv \
            --threads {threads} \
          > {log} 2>&1
        """
