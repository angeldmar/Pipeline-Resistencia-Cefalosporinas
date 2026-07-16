"""Pruebas unitarias de parse_quast.py (parte 11: evaluacion del ensamblaje)."""

import pandas as pd
import pytest

from parse_quast import classify_assembly, extract_quast_metrics, load_quast_report

DEFAULT_THRESHOLDS = dict(
    maximum_contigs=500, minimum_total_length=4_000_000, maximum_total_length=6_500_000,
    n50_warning_threshold=20_000,
)


def test_classify_assembly_pass():
    status = classify_assembly(contigs=80, total_length=5_000_000, n50=95_000, **DEFAULT_THRESHOLDS)
    assert status == "PASS"


def test_classify_assembly_warning_on_low_n50():
    status = classify_assembly(contigs=180, total_length=4_900_000, n50=12_000, **DEFAULT_THRESHOLDS)
    assert status == "WARNING"


@pytest.mark.parametrize(
    "contigs,total_length,n50",
    [
        (650, 5_100_000, 25_000),  # demasiados contigs, aunque N50 sea razonable
        (80, 3_000_000, 95_000),  # longitud total demasiado corta
        (80, 8_000_000, 95_000),  # longitud total demasiado larga
    ],
)
def test_classify_assembly_fail_takes_priority_over_warning(contigs, total_length, n50):
    # FAIL debe imponerse sobre WARNING incluso cuando el N50 es bueno: un
    # ensamblaje con demasiados contigs o longitud fuera de rango no es
    # confiable sin importar que tan continuo se vea localmente.
    status = classify_assembly(contigs=contigs, total_length=total_length, n50=n50, **DEFAULT_THRESHOLDS)
    assert status == "FAIL"


def build_quast_report_tsv(tmp_path, contigs=80, largest_contig=280_000, total_length=5_050_000, gc=50.65, n50=95_000):
    report_path = tmp_path / "report.tsv"
    report_path.write_text(
        "Assembly\tcontigs_filtered\n"
        "# contigs (>= 0 bp)\t120\n"
        f"# contigs\t{contigs}\n"
        f"Largest contig\t{largest_contig}\n"
        f"Total length\t{total_length}\n"
        f"GC (%)\t{gc}\n"
        f"N50\t{n50}\n"
    )
    return report_path


def test_load_quast_report_uses_exact_row_names_not_qualified_variants(tmp_path):
    # QUAST repite metricas calificadas por longitud (ej. "# contigs (>= 0 bp)")
    # ademas de la fila sin calificar ("# contigs"). extract_quast_metrics
    # debe usar la version sin calificar, no la primera coincidencia parcial.
    report_path = build_quast_report_tsv(tmp_path, contigs=80)
    report = load_quast_report(report_path)

    metrics = extract_quast_metrics("EC001", report, **DEFAULT_THRESHOLDS)

    assert metrics["contigs"] == 80  # no 120 (la fila calificada "(>= 0 bp)")


def test_extract_quast_metrics_end_to_end(tmp_path):
    report_path = build_quast_report_tsv(tmp_path)
    report = load_quast_report(report_path)

    metrics = extract_quast_metrics("EC001", report, **DEFAULT_THRESHOLDS)

    assert metrics == {
        "sample_id": "EC001",
        "contigs": 80,
        "largest_contig": 280_000,
        "total_length": 5_050_000,
        "gc_content_percent": 50.65,
        "n50": 95_000,
        "assembly_status": "PASS",
    }
