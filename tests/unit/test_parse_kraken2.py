"""Pruebas unitarias de parse_kraken2.py (parte 13: identificacion taxonomica)."""

import pandas as pd
import pytest

from parse_kraken2 import classify_taxonomy, extract_taxonomy_metrics

DEFAULT_THRESHOLDS = dict(minimum_ecoli_percentage=90, warning_ecoli_percentage=70, maximum_contaminant_percentage=5)


@pytest.mark.parametrize(
    "ecoli_pct,other_pct,expected",
    [
        (95, 3, "PASS"),
        (90, 4.99, "PASS"),
        (85, 3, "WARNING"),  # no llega al 90% de PASS pero si al 70% de WARNING
        (70, 3, "WARNING"),
        (69.9, 0, "FAIL"),
        (40, 45, "FAIL"),
    ],
)
def test_classify_taxonomy_thresholds(ecoli_pct, other_pct, expected):
    assert classify_taxonomy(ecoli_pct, other_pct, **DEFAULT_THRESHOLDS) == expected


def build_kraken2_report(rows: list[tuple[float, str, str]]) -> pd.DataFrame:
    """rows: lista de (percentage, rank_code, name)."""
    return pd.DataFrame(
        [{"percentage": pct, "reads_in_clade": 0, "reads_direct": 0, "rank_code": rank, "taxid": 0, "name": name}
         for pct, rank, name in rows]
    )


def test_shigella_excluded_from_contamination_but_flags_manual_review():
    # Caso central de la parte 13: una muestra con 85% E. coli + 10%
    # Shigella no debe tratar ese 10% como contaminacion (los cuenta aparte),
    # pero SI debe marcarse para revision manual.
    report = build_kraken2_report([
        (85.0, "S", "Escherichia coli"),
        (10.0, "S", "Shigella flexneri"),
        (3.0, "S", "Klebsiella pneumoniae"),
    ])

    metrics = extract_taxonomy_metrics("EC_SHIGELLA", report, **DEFAULT_THRESHOLDS)

    assert metrics["shigella_percentage"] == 10.0
    assert metrics["other_contaminant_percentage"] == 3.0  # Klebsiella solamente, sin Shigella
    assert metrics["requires_manual_review"] is True
    # 85% no alcanza el 90% de PASS, independientemente de Shigella.
    assert metrics["taxonomy_status"] == "WARNING"


def test_predominant_taxon_is_not_assumed_to_be_ecoli():
    report = build_kraken2_report([
        (40.0, "S", "Escherichia coli"),
        (45.0, "S", "Klebsiella pneumoniae"),
    ])

    metrics = extract_taxonomy_metrics("EC_FAIL", report, **DEFAULT_THRESHOLDS)

    assert metrics["predominant_taxon"] == "Klebsiella pneumoniae"
    assert metrics["taxonomy_status"] == "FAIL"


def test_no_species_level_rows_handled_without_crashing():
    # Caso negativo: un reporte sin ninguna fila a nivel de especie (ej. una
    # muestra totalmente sin clasificar) no debe hacer fallar el script.
    report = build_kraken2_report([(100.0, "U", "unclassified")])

    metrics = extract_taxonomy_metrics("EC_UNCLASSIFIED", report, **DEFAULT_THRESHOLDS)

    assert metrics["predominant_taxon"] == "unclassified"
    assert metrics["ecoli_percentage"] == 0.0
    assert metrics["taxonomy_status"] == "FAIL"
