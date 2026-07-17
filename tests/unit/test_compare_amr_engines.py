"""Pruebas unitarias de compare_amr_engines.py (parte 25: concordancia entre
motores de deteccion de AMR)."""

import pytest

from compare_amr_engines import (
    build_engine_agreement_long_table,
    build_engine_concordance_report,
    exact_gene_concordance,
    jaccard_similarity,
)


def test_exact_and_jaccard_concordance_functions():
    assert exact_gene_concordance({"blaCTX-M"}, {"blaCTX-M"}) == 1.0
    assert exact_gene_concordance({"blaCTX-M"}, {"blaTEM"}) == 0.0
    assert jaccard_similarity(set(), set()) == 1.0
    assert jaccard_similarity({"a", "b"}, {"a"}) == pytest.approx(0.5)


def test_concordance_report_detects_real_disagreement():
    # EC001: ambos motores concuerdan. EC002: solo AMRFinderPlus detecta
    # algo (discordancia real). EC003: ninguno detecta nada (concordante).
    amrfinder_families = {"EC001": {"blaCTX-M", "blaTEM"}, "EC002": {"blaCMY"}}
    abricate_families = {"EC001": {"blaCTX-M", "blaTEM"}}

    report = build_engine_concordance_report(
        ["EC001", "EC002", "EC003"], amrfinder_families, abricate_families
    ).set_index("sample_id")

    assert report.loc["EC001", "exact_gene_family_concordance"] == 1.0
    assert report.loc["EC001", "jaccard_similarity"] == 1.0

    assert report.loc["EC002", "exact_gene_family_concordance"] == 0.0
    assert report.loc["EC002", "jaccard_similarity"] == 0.0
    assert report.loc["EC002", "abricate_gene_families"] == "none"

    assert report.loc["EC003", "exact_gene_family_concordance"] == 1.0
    assert report.loc["EC003", "amrfinder_gene_families"] == "none"
    assert report.loc["EC003", "abricate_gene_families"] == "none"


def test_agreement_long_table_only_includes_families_detected_by_at_least_one_engine():
    amrfinder_families = {"EC001": {"blaCTX-M"}}
    abricate_families = {"EC001": {"blaCTX-M", "blaTEM"}}

    agreement_table = build_engine_agreement_long_table(["EC001"], amrfinder_families, abricate_families)

    assert len(agreement_table) == 2
    ctx_m_row = agreement_table.loc[agreement_table["gene_family"] == "blaCTX-M"].iloc[0]
    tem_row = agreement_table.loc[agreement_table["gene_family"] == "blaTEM"].iloc[0]

    assert ctx_m_row["amrfinder_result"] == "detected"
    assert ctx_m_row["abricate_result"] == "detected"
    assert tem_row["amrfinder_result"] == "not_detected"
    assert tem_row["abricate_result"] == "detected"


def test_agreement_long_table_empty_when_no_engine_detects_anything():
    agreement_table = build_engine_agreement_long_table(["EC003"], {}, {})
    assert agreement_table.empty
