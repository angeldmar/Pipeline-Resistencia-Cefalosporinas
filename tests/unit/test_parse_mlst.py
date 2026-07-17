"""Pruebas unitarias de parse_mlst.py (parte 26: tipificacion MLST)."""

import pandas as pd
import pytest

from parse_mlst import classify_sequence_type_status, normalize_mlst_output


@pytest.mark.parametrize(
    "sequence_type,allele_calls,expected_status",
    [
        ("131", ["adk(6)", "fumC(4)", "gyrB(4)"], "exact"),
        ("10", ["adk(~6)", "fumC(4)", "gyrB(4)"], "novel_allele"),  # coincidencia aproximada
        ("10", ["adk(6)", "fumC(?4)", "gyrB(4)"], "novel_allele"),
        ("-", ["adk(99)", "fumC(4)", "gyrB(4)"], "unresolved"),  # sin ST catalogado
    ],
)
def test_classify_sequence_type_status(sequence_type, allele_calls, expected_status):
    assert classify_sequence_type_status(sequence_type, allele_calls) == expected_status


def build_raw_mlst_table(scheme="ecoli", sequence_type="131", allele_calls=None) -> pd.DataFrame:
    allele_calls = allele_calls or ["adk(6)", "fumC(4)", "gyrB(4)", "icd(8)", "mdh(5)", "purA(3)", "recA(2)"]
    row = ["contigs.fasta", scheme, sequence_type] + allele_calls
    return pd.DataFrame([row])


def test_normalize_mlst_output_exact_match():
    raw_table = build_raw_mlst_table(sequence_type="131")

    result = normalize_mlst_output("EC001", raw_table)

    assert result["sample_id"] == "EC001"
    assert result["scheme"] == "ecoli"
    assert result["sequence_type"] == "131"
    assert result["sequence_type_status"] == "exact"
    assert "adk(6)" in result["allele_profile"]


def test_normalize_mlst_output_unresolved_sequence_type():
    raw_table = build_raw_mlst_table(sequence_type="-")

    result = normalize_mlst_output("EC_UNRESOLVED", raw_table)

    assert result["sequence_type"] == "-"
    assert result["sequence_type_status"] == "unresolved"


def test_normalize_mlst_output_handles_empty_table_without_crashing():
    # Caso negativo: si por alguna razon no hay fila de salida, no debe fallar
    # con un IndexError -- debe quedar marcado como no resuelto, visible.
    empty_table = pd.DataFrame()

    result = normalize_mlst_output("EC_EMPTY", empty_table)

    assert result["sequence_type_status"] == "unresolved"
    assert result["allele_profile"] == "none"
