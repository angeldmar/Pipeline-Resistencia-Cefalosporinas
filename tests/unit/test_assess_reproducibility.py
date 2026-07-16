"""Pruebas unitarias de assess_reproducibility.py (parte 19)."""

import pytest

from assess_reproducibility import (
    exact_gene_concordance,
    group_replicate_runs,
    jaccard_similarity,
    parse_replicate_run_id,
)


def test_exact_gene_concordance_matches_design_document_function():
    assert exact_gene_concordance({"blaCTX-M-15"}, {"blaCTX-M-15"}) == 1.0
    assert exact_gene_concordance({"blaCTX-M-15"}, {"blaCTX-M-27"}) == 0.0
    assert exact_gene_concordance(set(), set()) == 1.0


def test_jaccard_similarity_matches_design_document_function():
    assert jaccard_similarity(set(), set()) == 1.0
    assert jaccard_similarity({"a", "b"}, {"a", "b", "c"}) == pytest.approx(2 / 3)
    assert jaccard_similarity({"a"}, {"b"}) == 0.0
    assert jaccard_similarity({"a", "b"}, {"a", "b"}) == 1.0


@pytest.mark.parametrize(
    "sample_id,expected",
    [
        ("EC001_run1", ("EC001", 1)),
        ("EC001_run23", ("EC001", 23)),
        ("SAMPLE-WITH-DASHES_run2", ("SAMPLE-WITH-DASHES", 2)),
        ("EC001", None),  # muestra normal, sin sufijo de corrida
        ("EC001_runX", None),  # sufijo no numerico, no coincide con la convencion
    ],
)
def test_parse_replicate_run_id(sample_id, expected):
    assert parse_replicate_run_id(sample_id) == expected


def test_group_replicate_runs_orders_by_run_number_and_ignores_non_replicate_samples():
    sample_ids = ["EC001_run3", "EC001_run1", "EC001_run2", "EC002_run1", "EC003"]

    groups = group_replicate_runs(sample_ids)

    assert groups == {
        "EC001": ["EC001_run1", "EC001_run2", "EC001_run3"],
        "EC002": ["EC002_run1"],
    }
    assert "EC003" not in groups
