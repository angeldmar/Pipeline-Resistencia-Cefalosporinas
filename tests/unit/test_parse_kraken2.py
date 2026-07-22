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
    # pero SI debe marcarse para revision manual (muy por encima del umbral
    # de ruido por defecto, 0.1%).
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


def test_shigella_trace_below_noise_threshold_does_not_flag_manual_review():
    # Encontrado con datos reales (ERR17582235): Kraken2 asigna un rastro
    # minimo de lecturas a Shigella por ambiguedad de k-mers entre generos
    # cercanos, sin que sea una senal real de mezcla de especies. Por debajo
    # del umbral configurado (0.1% por defecto), NO debe disparar revision
    # manual -- si cualquier valor mayor a cero la disparara, la alerta
    # perderia utilidad al activarse casi siempre en datos reales.
    report = build_kraken2_report([
        (94.85, "S", "Escherichia coli"),
        (0.04, "S", "Shigella sonnei"),
    ])

    metrics = extract_taxonomy_metrics("EC_TRACE_SHIGELLA", report, **DEFAULT_THRESHOLDS)

    assert metrics["shigella_percentage"] == 0.04
    assert metrics["requires_manual_review"] is False


def test_shigella_review_threshold_is_configurable():
    report = build_kraken2_report([
        (94.0, "S", "Escherichia coli"),
        (0.5, "S", "Shigella sonnei"),
    ])

    # Con el umbral por defecto (0.1%), 0.5% si dispara revision.
    default_metrics = extract_taxonomy_metrics("EC_CONFIG_DEFAULT", report, **DEFAULT_THRESHOLDS)
    assert default_metrics["requires_manual_review"] is True

    # Con un umbral mas alto que la senal observada, no deberia dispararse.
    strict_metrics = extract_taxonomy_metrics(
        "EC_CONFIG_STRICT", report, **DEFAULT_THRESHOLDS, shigella_review_threshold_percentage=1.0
    )
    assert strict_metrics["requires_manual_review"] is False


def test_family_percentage_captures_unresolved_species_reads():
    # Encontrado con datos reales (ERR17582235, base Kraken2 recortada por
    # tamano): solo el 4.71% de las lecturas resolvio hasta especie, pero el
    # 90.52% resolvio al menos hasta familia Enterobacteriaceae (incluye el
    # 4.71% de especie mas lecturas que no bajaron mas alla de familia).
    # family_percentage debe capturar ese acumulado por separado de
    # ecoli_percentage -- se probo primero con nivel de genero, pero en la
    # practica casi nada se detiene justo ahi (ver docstring del modulo).
    report = build_kraken2_report([
        (90.52, "F", "Enterobacteriaceae"),
        (5.47, "G", "Escherichia"),
        (4.71, "S", "Escherichia coli"),
    ])

    metrics = extract_taxonomy_metrics("EC_FAMILY_ONLY", report, **DEFAULT_THRESHOLDS)

    assert metrics["ecoli_percentage"] == 4.71
    assert metrics["family_percentage"] == 90.52


def test_family_percentage_defaults_to_zero_without_family_level_row():
    report = build_kraken2_report([(95.0, "S", "Escherichia coli")])

    metrics = extract_taxonomy_metrics("EC_NO_FAMILY_ROW", report, **DEFAULT_THRESHOLDS)

    assert metrics["family_percentage"] == 0.0


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
