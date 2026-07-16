"""Pruebas negativas (seccion 22 del diseno del pipeline).

Varios casos negativos individuales (identificador duplicado, fenotipo
invalido, FASTQ/ensamblaje vacio, archivo corrupto) ya quedan cubiertos por
las pruebas unitarias de cada script (ver tests/unit/). Este archivo se
concentra en casos negativos que solo se manifiestan al CRUZAR modulos: una
muestra de "especie distinta" o "contaminada" no debe desaparecer del
pipeline, debe quedar EXCLUDED en la tabla maestra, visible y con motivo
identificable -- y algunos casos de entrada malformada que no se probaron
en otro lado.
"""

from pathlib import Path

import pandas as pd
import pytest

from generate_report import load_master_row
from merge_results import build_master_table
from validate_samples import validate_samples


def test_different_species_sample_is_excluded_not_silently_dropped():
    # "Especie distinta": Kraken2 clasifica la muestra como otra cosa
    # (taxonomy_status FAIL). No debe desaparecer de la tabla maestra.
    samples_table = pd.DataFrame([{"sample_id": "EC_WRONG_SPECIES", "expected_genes": "none"}])
    taxonomy_table = pd.DataFrame([{"sample_id": "EC_WRONG_SPECIES", "taxonomy_status": "FAIL"}])

    master_table = build_master_table(
        samples_table, fastp_table=None, quast_table=None, checkm_table=None,
        taxonomy_table=taxonomy_table, amr_table=None, reference_comparison_table=None,
    )

    assert len(master_table) == 1  # sigue presente
    assert master_table.iloc[0]["final_status"] == "EXCLUDED"
    assert master_table.iloc[0]["taxonomy_status"] == "FAIL"  # motivo identificable


def test_contaminated_sample_is_excluded_not_silently_dropped():
    # "Muestra contaminada": CheckM reporta contaminacion por encima del
    # umbral (completeness_status FAIL). Tampoco debe desaparecer.
    samples_table = pd.DataFrame([{"sample_id": "EC_CONTAMINATED", "expected_genes": "none"}])
    checkm_table = pd.DataFrame([{"sample_id": "EC_CONTAMINATED", "completeness_status": "FAIL"}])

    master_table = build_master_table(
        samples_table, fastp_table=None, quast_table=None, checkm_table=checkm_table,
        taxonomy_table=None, amr_table=None, reference_comparison_table=None,
    )

    assert len(master_table) == 1
    assert master_table.iloc[0]["final_status"] == "EXCLUDED"


def test_missing_sample_in_master_table_raises_clear_error(tmp_path):
    # Pedir el reporte de una muestra que no existe en master_results.tsv no
    # debe fallar con un error críptico de pandas, sino con un mensaje claro
    # que identifique la muestra faltante.
    master_results_path = tmp_path / "master_results.tsv"
    pd.DataFrame([{"sample_id": "EC001", "final_status": "PASS"}]).to_csv(
        master_results_path, sep="\t", index=False
    )

    with pytest.raises(ValueError, match="EC999"):
        load_master_row(master_results_path, "EC999")


def test_samples_file_with_zero_rows_is_accepted_but_produces_no_samples(tmp_path):
    # Metadatos "incompletos" en su forma mas extrema: un samples.tsv con
    # encabezado pero cero muestras. No debe fallar la validacion (las
    # columnas requeridas si estan), pero tampoco debe inventar muestras.
    samples_path = tmp_path / "samples.tsv"
    header = "\t".join([
        "sample_id", "run_accession", "biosample", "sequencing_platform",
        "phenotype_cefotaxime", "phenotype_ceftriaxone", "phenotype_ceftazidime",
        "expected_genes", "data_source",
    ])
    samples_path.write_text(header + "\n")

    validated_table = validate_samples(samples_path)

    assert len(validated_table) == 0


def test_samples_file_missing_multiple_required_columns_reports_all_of_them(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    # Solo quedan 3 de las 9 columnas requeridas.
    samples_path.write_text("sample_id\trun_accession\tbiosample\nEC001\tSRR000001\tSAMN000001\n")

    with pytest.raises(ValueError) as exc_info:
        validate_samples(samples_path)

    error_message = str(exc_info.value)
    for missing_column in ["sequencing_platform", "phenotype_cefotaxime", "data_source"]:
        assert missing_column in error_message
