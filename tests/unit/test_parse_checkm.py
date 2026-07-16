"""Pruebas unitarias de parse_checkm.py (parte 12: completitud y contaminacion)."""

import pytest

from parse_checkm import build_exclusion_registry, classify_completeness
import pandas as pd


@pytest.mark.parametrize(
    "completeness,contamination,expected_status",
    [
        (98.7, 1.2, "PASS"),
        (95.0, 4.99, "PASS"),
        (94.9, 1.0, "FAIL"),  # completitud justo debajo del minimo
        (99.0, 5.0, "FAIL"),  # contaminacion justo en el maximo (no permitido)
        (89.5, 0.5, "FAIL"),  # solo falla por completitud
        (99.0, 6.3, "FAIL"),  # solo falla por contaminacion
    ],
)
def test_classify_completeness(completeness, contamination, expected_status):
    status = classify_completeness(completeness, contamination, minimum_completeness=95, maximum_contamination=5)
    assert status == expected_status


def test_exclusion_registry_only_lists_failed_samples_with_reason():
    combined_table = pd.DataFrame([
        {"sample_id": "EC_PASS", "completeness_percent": 98.7, "contamination_percent": 1.2, "completeness_status": "PASS"},
        {"sample_id": "EC_FAIL", "completeness_percent": 89.5, "contamination_percent": 0.5, "completeness_status": "FAIL"},
    ])

    exclusions = build_exclusion_registry(combined_table)

    assert list(exclusions["sample_id"]) == ["EC_FAIL"]
    assert (exclusions["reason"] == "completitud/contaminacion fuera de umbral (CheckM)").all()
