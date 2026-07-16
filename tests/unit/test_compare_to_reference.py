"""Pruebas unitarias de compare_to_reference.py (parte 16).

Incluye el caso que corrigio un bug real durante el desarrollo: una
referencia indeterminada (None) no debe caer silenciosamente en la rama de
"referencia negativa" solo porque None tambien es falsy en Python.
"""

import pandas as pd
import pytest

from compare_to_reference import (
    build_comparison_row,
    confusion_category,
    determine_reference_status,
)


@pytest.mark.parametrize(
    "reference_positive,pipeline_positive,expected",
    [
        (True, True, "TP"),
        (False, False, "TN"),
        (False, True, "FP"),
        (True, False, "FN"),
    ],
)
def test_confusion_category_matches_design_document_function(reference_positive, pipeline_positive, expected):
    assert confusion_category(reference_positive, pipeline_positive) == expected


def test_confusion_category_indeterminate_when_reference_is_none():
    # Quinta categoria (mas alla del ejemplo del documento): sin estandar de
    # referencia documentado, no se puede clasificar como TP/TN/FP/FN.
    assert confusion_category(None, True) == "Indeterminado"
    assert confusion_category(None, False) == "Indeterminado"


@pytest.mark.parametrize(
    "expected_gene,expected_positive,expected_status",
    [
        ("blaCTX-M-15", True, "positive"),
        ("none", False, "negative"),
        ("NA", None, "indeterminate"),
        ("", None, "indeterminate"),
        ("nan", None, "indeterminate"),
    ],
)
def test_determine_reference_status(expected_gene, expected_positive, expected_status):
    reference_positive, reference_status = determine_reference_status(expected_gene)
    assert reference_positive == expected_positive
    assert reference_status == expected_status


def make_amr_rows(rows: list[dict]) -> pd.DataFrame:
    # Se fijan las columnas explicitamente para que una lista vacia (una
    # muestra sin genes detectados) siga produciendo un DataFrame con las
    # columnas esperadas, en vez de uno completamente vacio sin columnas.
    return pd.DataFrame(rows, columns=["gene_symbol", "antimicrobial_class"])


def test_exact_match_is_true_positive():
    sample_rows = make_amr_rows([
        {"gene_symbol": "blaCTX-M-15", "antimicrobial_class": "BETA-LACTAM"},
    ])
    result = build_comparison_row("EC001", "blaCTX-M-15", sample_rows)
    assert result["match_type"] == "exact"
    assert result["confusion_category"] == "TP"


def test_family_match_with_different_allele_still_counts_as_detected():
    # Se esperaba blaCTX-M-15 pero se detecto blaCTX-M-27 (misma familia).
    sample_rows = make_amr_rows([
        {"gene_symbol": "blaCTX-M-27", "antimicrobial_class": "BETA-LACTAM"},
    ])
    result = build_comparison_row("EC001", "blaCTX-M-15", sample_rows)
    assert result["match_type"] == "family"
    assert result["detected_gene"] == "blaCTX-M-27"
    assert result["confusion_category"] == "TP"


def test_expected_gene_not_detected_is_false_negative():
    sample_rows = make_amr_rows([])
    result = build_comparison_row("EC004", "blaCTX-M-27", sample_rows)
    assert result["confusion_category"] == "FN"
    assert result["detected_gene"] == "none"


def test_no_gene_expected_and_none_detected_is_true_negative():
    sample_rows = make_amr_rows([])
    result = build_comparison_row("EC002", "none", sample_rows)
    assert result["confusion_category"] == "TN"


def test_unexpected_beta_lactam_gene_is_false_positive():
    sample_rows = make_amr_rows([
        {"gene_symbol": "blaTEM-1", "antimicrobial_class": "BETA-LACTAM"},
    ])
    result = build_comparison_row("EC005", "none", sample_rows)
    assert result["confusion_category"] == "FP"
    assert result["detected_gene"] == "blaTEM-1"


def test_indeterminate_reference_reports_detected_genes_without_scoring():
    # Regresion del bug encontrado durante las pruebas manuales: antes,
    # reference_positive=None caia en la rama de "no se espera nada" (por
    # ser falsy) y calificaba genes como "unexpected". Ahora debe usar la
    # rama dedicada, con match_type="not_applicable".
    sample_rows = make_amr_rows([
        {"gene_symbol": "blaCTX-M-27", "antimicrobial_class": "BETA-LACTAM"},
    ])
    result = build_comparison_row("EC006", "NA", sample_rows)

    assert result["reference_status"] == "indeterminate"
    assert result["confusion_category"] == "Indeterminado"
    assert result["match_type"] == "not_applicable"
    assert result["detected_gene"] == "blaCTX-M-27"  # se informa igual, sin descartarlo


def test_indeterminate_reference_with_no_detections_reports_none():
    sample_rows = make_amr_rows([])
    result = build_comparison_row("EC007", "NA", sample_rows)
    assert result["confusion_category"] == "Indeterminado"
    assert result["detected_gene"] == "none"
