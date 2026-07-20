"""Pruebas unitarias de parse_fastp.py (partes 8 y 9: calidad y cobertura)."""

import pytest

from parse_fastp import classify_coverage, estimate_coverage, extract_fastp_metrics


def test_estimate_coverage_matches_reference_example():
    # Caso de referencia: 2,000,000
    # lecturas x 150 pb / 5,000,000 pb de genoma = 60x.
    result = estimate_coverage(read_count=2_000_000, mean_read_length=150, genome_size=5_000_000)
    assert result == 60


def test_estimate_coverage_rejects_non_positive_genome_size():
    with pytest.raises(ValueError):
        estimate_coverage(read_count=1000, mean_read_length=150, genome_size=0)


@pytest.mark.parametrize(
    "coverage,expected_status",
    [
        (60, "PASS"),
        (30, "PASS"),
        (29.9, "WARNING"),
        (15, "WARNING"),
        (14.9, "FAIL"),
        (0, "FAIL"),
    ],
)
def test_classify_coverage_thresholds(coverage, expected_status):
    assert classify_coverage(coverage, minimum_coverage=30, warning_coverage=15) == expected_status


def build_fastp_report(initial_reads: int, retained_reads: int, duplication_rate: float, gc_content: float) -> dict:
    total_bases = retained_reads * 150
    return {
        "summary": {
            "before_filtering": {"total_reads": initial_reads, "total_bases": initial_reads * 150},
            "after_filtering": {
                "total_reads": retained_reads,
                "total_bases": total_bases,
                "q20_bases": int(total_bases * 0.98),
                "q30_bases": int(total_bases * 0.95),
                "gc_content": gc_content,
            },
        },
        "duplication": {"rate": duplication_rate},
    }


def test_extract_fastp_metrics_computes_coverage_and_status():
    report = build_fastp_report(initial_reads=2_001_000, retained_reads=2_000_000, duplication_rate=0.02, gc_content=0.50)

    metrics = extract_fastp_metrics(
        "EC001", report, genome_size=5_000_000, minimum_coverage=30, warning_coverage=15
    )

    assert metrics["sample_id"] == "EC001"
    assert metrics["initial_reads"] == 2_001_000
    assert metrics["retained_reads"] == 2_000_000
    assert metrics["mean_read_length"] == 150.0
    assert metrics["estimated_coverage"] == 60.0
    assert metrics["coverage_status"] == "PASS"
    assert metrics["gc_content_percent"] == 50.0
    assert metrics["duplication_rate_percent"] == 2.0


def test_extract_fastp_metrics_handles_zero_retained_reads_without_crashing():
    # Caso negativo: un FASTQ vacio (o completamente filtrado) no debe
    # producir una division por cero ni un traceback -- debe marcarse FAIL.
    report = build_fastp_report(initial_reads=1000, retained_reads=0, duplication_rate=0.0, gc_content=0.0)

    metrics = extract_fastp_metrics(
        "EC_EMPTY", report, genome_size=5_000_000, minimum_coverage=30, warning_coverage=15
    )

    assert metrics["retained_reads"] == 0
    assert metrics["estimated_coverage"] == 0.0
    assert metrics["coverage_status"] == "FAIL"
