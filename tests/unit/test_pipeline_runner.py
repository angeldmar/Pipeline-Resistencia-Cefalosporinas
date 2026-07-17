"""Pruebas unitarias de webapp/pipeline_runner.py (parte 28: interfaz web
local de analisis ad-hoc)."""

import pytest

import pipeline_runner as runner


@pytest.mark.parametrize(
    "sample_id",
    ["EC001", "EC_TEST-01", "a", "A" * 64],
)
def test_validate_sample_id_accepts_safe_identifiers(sample_id):
    runner.validate_sample_id(sample_id)  # no debe lanzar


@pytest.mark.parametrize(
    "sample_id",
    [
        "",  # vacio
        "../../etc/passwd",  # traversal de ruta
        "EC 001",  # espacio
        "EC/001",  # separador de ruta
        "A" * 65,  # demasiado largo
        "EC;rm -rf /",  # intento de inyeccion de shell
    ],
)
def test_validate_sample_id_rejects_unsafe_identifiers(sample_id):
    with pytest.raises(runner.InvalidSampleIdError):
        runner.validate_sample_id(sample_id)


def test_create_isolated_samples_file_never_touches_main_samples_tsv(tmp_path, monkeypatch, repo_root):
    # La tabla aislada debe ser un archivo propio de esta carga, y nunca
    # modificar config/samples.tsv (el lote curado principal).
    monkeypatch.setattr(runner, "RUNS_DIR", tmp_path)

    samples_path = runner.create_isolated_samples_file("EC_TEST_SAMPLE", "ILLUMINA", "")

    assert samples_path != repo_root / "config" / "samples.tsv"
    assert samples_path.is_file()
    assert samples_path.parent == tmp_path / "EC_TEST_SAMPLE"


def test_create_isolated_samples_file_defaults_expected_genes_to_na(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "RUNS_DIR", tmp_path)

    samples_path = runner.create_isolated_samples_file("EC_TEST_SAMPLE", "ILLUMINA", "")

    content = samples_path.read_text()
    assert "\tNA\tCarga local ad-hoc (interfaz web)\n" in content


def test_create_isolated_samples_file_preserves_provided_expected_gene(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "RUNS_DIR", tmp_path)

    samples_path = runner.create_isolated_samples_file("EC_TEST_SAMPLE", "ILLUMINA", "blaCTX-M-15")

    content = samples_path.read_text()
    assert "\tblaCTX-M-15\tCarga local ad-hoc (interfaz web)\n" in content


def test_status_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "RUNS_DIR", tmp_path)

    assert runner.read_status("EC_NEW_SAMPLE") == "not_found"

    runner.write_status("EC_NEW_SAMPLE", "running")
    assert runner.read_status("EC_NEW_SAMPLE") == "running"

    runner.write_status("EC_NEW_SAMPLE", "done")
    assert runner.read_status("EC_NEW_SAMPLE") == "done"
