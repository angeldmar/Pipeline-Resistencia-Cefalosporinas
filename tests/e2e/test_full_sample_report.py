"""Prueba de extremo a extremo: de los resultados de cada modulo hasta el
reporte HTML final de una muestra.

Alcance: procesa una muestra conocida y sintetica a traves de TODA la
cadena de scripts de este pipeline en Python (parseo -> clasificacion ->
comparacion -> integracion -> generacion de reporte), partiendo de
artefactos con el formato real que producirian las herramientas externas
(JSON de fastp, report.tsv de QUAST, tabla de CheckM, reporte de Kraken2,
tabla de AMRFinderPlus). No incluye correr las herramientas bioinformaticas
en si (no estan disponibles en este entorno de desarrollo, ver la parte de
ambientes Conda), pero SI ejercita cada linea de codigo Python real del
pipeline, con datos que imitan exactamente el formato de esas herramientas.
"""

from pathlib import Path

import pandas as pd
import yaml

from classify_cephalosporin_genes import classify_amr_table, load_resistance_targets
from compare_to_reference import build_comparison_table
from generate_report import (
    build_gene_interpretation_sentences,
    collect_warnings,
    render_report,
)
from merge_results import build_master_table
from parse_amrfinder import normalize_amrfinder_table
from parse_checkm import extract_checkm_metrics, load_checkm_report
from parse_fastp import extract_fastp_metrics
from parse_kraken2 import extract_taxonomy_metrics, load_kraken2_report
from parse_quast import extract_quast_metrics, load_quast_report


def test_full_pipeline_from_tool_outputs_to_html_report(tmp_path, repo_root):
    sample_id = "EC001"

    # --- 1. fastp: JSON con formato real -----------------------------------
    fastp_report = {
        "summary": {
            "before_filtering": {"total_reads": 2_001_000, "total_bases": 2_001_000 * 150},
            "after_filtering": {
                "total_reads": 2_000_000, "total_bases": 2_000_000 * 150,
                "q20_bases": int(2_000_000 * 150 * 0.98), "q30_bases": int(2_000_000 * 150 * 0.95),
                "gc_content": 0.501,
            },
        },
        "duplication": {"rate": 0.02},
    }
    fastp_metrics = extract_fastp_metrics(
        sample_id, fastp_report, genome_size=5_000_000, minimum_coverage=30, warning_coverage=15
    )
    assert fastp_metrics["coverage_status"] == "PASS"

    # --- 2. QUAST: report.tsv con formato real ------------------------------
    quast_report_path = tmp_path / "quast_report.tsv"
    quast_report_path.write_text(
        "Assembly\tcontigs_filtered\n"
        "# contigs (>= 0 bp)\t120\n"
        "# contigs\t80\n"
        "Largest contig\t280000\n"
        "Total length\t5050000\n"
        "GC (%)\t50.65\n"
        "N50\t95000\n"
    )
    quast_metrics = extract_quast_metrics(
        sample_id, load_quast_report(quast_report_path),
        maximum_contigs=500, minimum_total_length=4_000_000, maximum_total_length=6_500_000,
        n50_warning_threshold=20_000,
    )
    assert quast_metrics["assembly_status"] == "PASS"

    # --- 3. CheckM: tabla --tab_table con formato real ----------------------
    checkm_report_path = tmp_path / "checkm_summary.tsv"
    checkm_report_path.write_text(
        "Bin Id\tMarker lineage\tCompleteness\tContamination\n"
        f"{sample_id}\tk__Bacteria\t98.7\t1.2\n"
    )
    checkm_metrics = extract_checkm_metrics(
        sample_id, load_checkm_report(checkm_report_path),
        minimum_completeness=95, maximum_contamination=5,
    )
    assert checkm_metrics["completeness_status"] == "PASS"

    # --- 4. Kraken2: reporte jerarquico con formato real --------------------
    kraken2_report_path = tmp_path / "kraken2_report.tsv"
    kraken2_report_path.write_text(
        " 95.00\t950000\t100\tD\t2\tBacteria\n"
        " 95.00\t950000\t950000\tS\t562\tEscherichia coli\n"
        "  3.00\t30000\t30000\tS\t573\tKlebsiella pneumoniae\n"
    )
    taxonomy_metrics = extract_taxonomy_metrics(
        sample_id, load_kraken2_report(kraken2_report_path),
        minimum_ecoli_percentage=90, warning_ecoli_percentage=70, maximum_contaminant_percentage=5,
    )
    assert taxonomy_metrics["taxonomy_status"] == "PASS"

    # --- 5. AMRFinderPlus: tabla con formato real ---------------------------
    raw_amrfinder_table = pd.DataFrame([{
        "Element symbol": "blaCTX-M-15", "Element name": "CTX-M-15 family class A ESBL",
        "Class": "BETA-LACTAM", "Subclass": "CEPHALOSPORIN", "Method": "ALLELEX",
        "% Identity to reference": 100.0, "% Coverage of reference": 100.0,
        "Contig id": "contig_1", "Start": 1, "Stop": 900,
    }])
    normalized_amr_table = normalize_amrfinder_table(
        sample_id, raw_amrfinder_table, minimum_identity=90, minimum_gene_coverage=80
    )
    resistance_targets = load_resistance_targets(repo_root / "config" / "resistance_targets.yaml")
    classified_amr_table = classify_amr_table(normalized_amr_table, resistance_targets)
    assert classified_amr_table.iloc[0]["beta_lactamase_category"] == "ESBL"

    # --- 6. Comparacion con el estandar de referencia -----------------------
    samples_table = pd.DataFrame([{"sample_id": sample_id, "expected_genes": "blaCTX-M-15"}])
    reference_comparison_table = build_comparison_table(samples_table, normalized_amr_table)
    assert reference_comparison_table.iloc[0]["confusion_category"] == "TP"

    # --- 7. Integracion: tabla maestra --------------------------------------
    master_table = build_master_table(
        samples_table,
        fastp_table=pd.DataFrame([fastp_metrics]),
        quast_table=pd.DataFrame([quast_metrics]),
        checkm_table=pd.DataFrame([checkm_metrics]),
        taxonomy_table=pd.DataFrame([taxonomy_metrics]),
        amr_table=normalized_amr_table,
        reference_comparison_table=reference_comparison_table,
        performance_by_sample_table=pd.DataFrame([
            {"sample_id": sample_id, "total_elapsed_seconds": 682.6, "peak_max_ram_gb": 7.81}
        ]),
    )
    assert master_table.iloc[0]["final_status"] == "PASS"

    # --- 8. Reporte HTML final -----------------------------------------------
    master_row = master_table.iloc[0].to_dict()
    amr_genes = classified_amr_table.to_dict("records")

    html_content = render_report(
        {
            **master_row,
            "amr_genes": amr_genes,
            "gene_interpretations": build_gene_interpretation_sentences(amr_genes),
            "warnings": collect_warnings(master_row, amr_genes),
            "amr_chart_base64": None,
            "tool_versions": {"fastp": "1.0.0"},
            "performance_by_module": [],
            "generated_at": "2026-07-16 00:00 UTC",
        },
        repo_root / "workflow" / "templates",
        "sample_report.html.j2",
    )

    # El reporte debe reflejar el hallazgo genotipico real...
    assert "blaCTX-M-15" in html_content
    assert "beta-lactamasas de espectro extendido" in html_content
    # ...sin emitir ninguna conclusion clinica prohibida.
    forbidden_phrases = ["aislado resistente", "resistente a ceftriaxona", "resistente a cefotaxima"]
    assert not any(phrase in html_content.lower() for phrase in forbidden_phrases)
    assert "PASS" in html_content
