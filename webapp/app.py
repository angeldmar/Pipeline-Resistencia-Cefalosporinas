"""Interfaz web local de analisis ad-hoc: subir un FASTQ paired-end o un
FASTA ya ensamblado y ver el reporte que el pipeline genera para el.

Herramienta de conveniencia, DISTINTA del pipeline en si (ver docstring de
pipeline_runner.py para el detalle de por que existen dos modos de
orquestacion). Requiere que las herramientas bioinformaticas del pipeline
esten instaladas (via los ambientes Conda de workflow/envs/, o de otra
forma accesibles en el PATH) para que el analisis real corra -- esta
interfaz no reemplaza esa instalacion, solo evita tener que editar
config/samples.tsv y correr Snakemake a mano para una muestra suelta.

Uso:
    pip install -r webapp/requirements.txt
    python webapp/app.py
    (abrir http://127.0.0.1:5000 en el navegador)
"""

from __future__ import annotations

from flask import Flask, redirect, render_template, request, send_file, url_for

import pipeline_runner as runner

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 500 * 1024 * 1024  # 500 MB: limite razonable para un FASTQ/FASTA de un genoma bacteriano


@app.route("/")
def index():
    return render_template("upload_form.html")


@app.route("/upload", methods=["POST"])
def upload():
    sample_id = request.form.get("sample_id", "").strip()
    mode = request.form.get("mode", "")
    sequencing_platform = request.form.get("sequencing_platform", "ILLUMINA").strip()
    expected_genes = request.form.get("expected_genes", "").strip()
    try:
        threads = max(1, int(request.form.get("threads", "4")))
    except ValueError:
        threads = 4

    try:
        runner.validate_sample_id(sample_id)
    except runner.InvalidSampleIdError as error:
        return render_template("upload_form.html", error=str(error)), 400

    if runner.report_path(sample_id).is_file() or runner.read_status(sample_id) == "running":
        return render_template(
            "upload_form.html",
            error=f"Ya existe una carga o un resultado para '{sample_id}'. Usa otro identificador.",
        ), 400

    if mode == "fastq":
        r1 = request.files.get("fastq_r1")
        r2 = request.files.get("fastq_r2")
        if not r1 or not r2 or not r1.filename or not r2.filename:
            return render_template("upload_form.html", error="Se necesitan ambos archivos FASTQ (R1 y R2)."), 400
        runner.save_uploaded_fastq(sample_id, r1, r2)
        runner.launch_fastq_pipeline(sample_id, sequencing_platform, expected_genes, threads)

    elif mode == "fasta":
        fasta = request.files.get("fasta_file")
        if not fasta or not fasta.filename:
            return render_template("upload_form.html", error="Se necesita un archivo FASTA."), 400
        assembly_path = runner.save_uploaded_fasta(sample_id, fasta)
        runner.run_fasta_only_pipeline(sample_id, assembly_path, expected_genes, threads)

    else:
        return render_template("upload_form.html", error="Selecciona un modo de carga (FASTQ o FASTA)."), 400

    return redirect(url_for("status", sample_id=sample_id))


@app.route("/status/<sample_id>")
def status(sample_id: str):
    try:
        runner.validate_sample_id(sample_id)
    except runner.InvalidSampleIdError as error:
        return str(error), 400

    current_status = runner.read_status(sample_id)
    log_path = runner.log_file(sample_id)
    log_tail = ""
    if log_path.is_file():
        log_lines = log_path.read_text(errors="replace").splitlines()
        log_tail = "\n".join(log_lines[-60:])  # solo las ultimas lineas, el log completo puede ser largo

    return render_template(
        "status.html",
        sample_id=sample_id,
        status=current_status,
        log_tail=log_tail,
        report_ready=runner.report_path(sample_id).is_file(),
    )


@app.route("/reports/<sample_id>")
def view_report(sample_id: str):
    try:
        runner.validate_sample_id(sample_id)
    except runner.InvalidSampleIdError as error:
        return str(error), 400

    path = runner.report_path(sample_id)
    if not path.is_file():
        return f"Todavia no hay reporte para '{sample_id}'.", 404
    return send_file(path)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
