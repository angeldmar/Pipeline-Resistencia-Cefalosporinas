"""Pruebas unitarias de generate_report.py (parte 27: descarga CSV desde el
reporte HTML, y partes 19/21: interpretacion genotipica, advertencias)."""

import base64

import pandas as pd

from generate_report import build_csv_download_data_uri, collect_warnings


def test_build_csv_download_data_uri_roundtrips_to_valid_csv():
    rows = [
        {"gene_symbol": "blaCTX-M-15", "percent_identity": 100.0},
        {"gene_symbol": "blaTEM-1", "percent_identity": 99.5},
    ]

    data_uri = build_csv_download_data_uri(rows)

    assert data_uri.startswith("data:text/csv;charset=utf-8;base64,")
    encoded_payload = data_uri.split(",", maxsplit=1)[1]
    decoded_csv = base64.b64decode(encoded_payload).decode("utf-8")

    # El CSV decodificado debe ser identico a lo que produciria pandas
    # directamente, confirmando que no se perdio ni corrompio informacion
    # en el viaje de ida y vuelta por base64.
    expected_csv = pd.DataFrame(rows).to_csv(index=False)
    assert decoded_csv == expected_csv


def test_build_csv_download_data_uri_returns_none_for_empty_rows():
    assert build_csv_download_data_uri([]) is None


def test_collect_warnings_flags_low_confidence_detections():
    master_row = {"final_status": "PASS"}
    amr_genes = [
        {"gene_symbol": "blaCTX-M-15", "meets_identity_coverage_threshold": True},
        {"gene_symbol": "aac(3)-IId", "meets_identity_coverage_threshold": False},
    ]

    warnings = collect_warnings(master_row, amr_genes)

    assert any("aac(3)-IId" in warning for warning in warnings)


def test_collect_warnings_explains_low_species_high_family_pattern():
    # Encontrado con datos reales (ERR17582235): taxonomy_status en WARNING/
    # FAIL por bajo % de especie, con un % de familia Enterobacteriaceae
    # bastante mas alto, sugiere limitacion de resolucion de la base de
    # datos, no una mezcla real de especies -- el reporte debe explicar
    # esto, no solo mostrar FAIL.
    master_row = {
        "final_status": "PASS", "taxonomy_status": "FAIL",
        "ecoli_percentage": 4.71, "family_percentage": 90.52,
    }

    warnings = collect_warnings(master_row, [])

    assert any("limitación de resolución" in warning for warning in warnings)


def test_collect_warnings_does_not_explain_pattern_when_gap_is_small():
    # Si el % de familia no es notablemente mas alto que el de especie, no
    # aplica la explicacion de resolucion de base de datos (la diferencia
    # es solo ruido normal, no un patron a destacar).
    master_row = {
        "final_status": "PASS", "taxonomy_status": "WARNING",
        "ecoli_percentage": 80.0, "family_percentage": 82.0,
    }

    warnings = collect_warnings(master_row, [])

    assert not any("limitación de resolución" in warning for warning in warnings)
