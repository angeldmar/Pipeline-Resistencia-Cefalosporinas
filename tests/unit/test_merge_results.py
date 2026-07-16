"""Pruebas unitarias de merge_results.py (parte 17: tabla maestra)."""

import pandas as pd
import pytest

from merge_results import build_amr_sample_summary, determine_final_status


@pytest.mark.parametrize(
    "statuses,expected",
    [
        (["PASS", "PASS", "PASS", "PASS"], "PASS"),
        (["PASS", "WARNING", "PASS", "PASS"], "WARNING"),
        (["PASS", "FAIL", "WARNING", "PASS"], "EXCLUDED"),  # FAIL manda, aunque haya WARNING tambien
        ([], "PENDING"),
        ([float("nan"), float("nan")], "PENDING"),
    ],
)
def test_determine_final_status(statuses, expected):
    assert determine_final_status(statuses) == expected


def test_amr_sample_summary_counts_only_confident_detections():
    amr_table = pd.DataFrame([
        {"sample_id": "EC001", "gene_symbol": "blaCTX-M-15", "antimicrobial_class": "BETA-LACTAM", "meets_identity_coverage_threshold": True},
        {"sample_id": "EC001", "gene_symbol": "aac(3)-IId", "antimicrobial_class": "AMINOGLYCOSIDE", "meets_identity_coverage_threshold": False},
        {"sample_id": "EC002", "gene_symbol": "blaCMY-2", "antimicrobial_class": "BETA-LACTAM", "meets_identity_coverage_threshold": True},
    ])

    summary = build_amr_sample_summary(amr_table).set_index("sample_id")

    assert summary.loc["EC001", "detected_gene_count"] == 1
    assert summary.loc["EC001", "detected_genes"] == "blaCTX-M-15"
    assert summary.loc["EC001", "detected_beta_lactam_gene_count"] == 1
    assert summary.loc["EC002", "detected_gene_count"] == 1
