"""Pruebas unitarias de validate_samples.py (parte 3)."""

import pandas as pd
import pytest

from validate_samples import (
    check_accession_formats,
    check_documented_source,
    check_duplicate_sample_ids,
    check_required_columns,
    check_valid_phenotypes,
    validate_samples,
)

VALID_ROW = {
    "sample_id": "EC001",
    "run_accession": "SRR000001",
    "biosample": "SAMN000001",
    "sequencing_platform": "ILLUMINA",
    "phenotype_cefotaxime": "R",
    "phenotype_ceftriaxone": "R",
    "phenotype_ceftazidime": "R",
    "expected_genes": "blaCTX-M-15",
    "data_source": "NCBI Pathogen Detection",
}


def make_samples_table(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


def test_valid_table_has_no_errors():
    table = make_samples_table([VALID_ROW])
    assert check_required_columns(table) == []
    assert check_duplicate_sample_ids(table) == []
    assert check_valid_phenotypes(table) == []
    assert check_accession_formats(table) == []
    assert check_documented_source(table) == []


def test_missing_required_column_detected():
    table = make_samples_table([VALID_ROW]).drop(columns=["data_source"])
    errors = check_required_columns(table)
    assert len(errors) == 1
    assert "data_source" in errors[0]


def test_duplicate_sample_id_detected():
    table = make_samples_table([VALID_ROW, {**VALID_ROW, "run_accession": "SRR000002"}])
    errors = check_duplicate_sample_ids(table)
    assert len(errors) == 1
    assert "EC001" in errors[0]


@pytest.mark.parametrize("invalid_phenotype", ["X", "Resistant", "1", ""])
def test_invalid_phenotype_detected(invalid_phenotype):
    row = {**VALID_ROW, "phenotype_cefotaxime": invalid_phenotype}
    table = make_samples_table([row])
    errors = check_valid_phenotypes(table)
    assert len(errors) == 1
    assert "phenotype_cefotaxime" in errors[0]


@pytest.mark.parametrize("valid_phenotype", ["S", "I", "R", "NA"])
def test_valid_phenotype_values_accepted(valid_phenotype):
    row = {**VALID_ROW, "phenotype_cefotaxime": valid_phenotype}
    table = make_samples_table([row])
    assert check_valid_phenotypes(table) == []


@pytest.mark.parametrize(
    "bad_run_accession", ["SRRBAD", "12345", "ERR", "srr000001"]
)
def test_malformed_run_accession_detected(bad_run_accession):
    row = {**VALID_ROW, "run_accession": bad_run_accession}
    table = make_samples_table([row])
    errors = check_accession_formats(table)
    assert len(errors) == 1
    assert "run_accession" in errors[0]


def test_malformed_biosample_detected():
    row = {**VALID_ROW, "biosample": "NOTABIOSAMPLE"}
    table = make_samples_table([row])
    errors = check_accession_formats(table)
    assert len(errors) == 1
    assert "biosample" in errors[0]


def test_missing_data_source_detected():
    row = {**VALID_ROW, "data_source": "NA"}
    table = make_samples_table([row])
    errors = check_documented_source(table)
    assert len(errors) == 1
    assert "EC001" in errors[0]


def test_validate_samples_raises_on_multiple_problems(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    bad_row = {**VALID_ROW, "phenotype_cefotaxime": "X", "data_source": "NA"}
    make_samples_table([bad_row]).to_csv(samples_path, sep="\t", index=False)

    with pytest.raises(ValueError) as exc_info:
        validate_samples(samples_path)

    # El error debe reportar AMBOS problemas de una sola vez, no solo el primero.
    assert "phenotype_cefotaxime" in str(exc_info.value)
    assert "fuente de datos" in str(exc_info.value)


def test_validate_samples_accepts_valid_file(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    make_samples_table([VALID_ROW]).to_csv(samples_path, sep="\t", index=False)

    validated_table = validate_samples(samples_path)

    assert list(validated_table["sample_id"]) == ["EC001"]
