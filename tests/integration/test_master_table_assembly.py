"""Prueba de integracion: ensamblaje de la tabla maestra a partir de TODOS
los modulos (parte 17), incluyendo el desempeno (parte 18).

Formaliza como pytest la prueba manual de extremo a extremo que se hizo
durante el desarrollo de merge_results.py: tres muestras consistentes que
pasan por fastp, QUAST, CheckM, Kraken2, AMRFinderPlus, comparacion de
referencia y desempeno, verificando que la tabla maestra final combine todo
correctamente.
"""

import pandas as pd

from merge_results import build_master_table


def test_master_table_combines_all_modules_and_computes_final_status():
    samples_table = pd.DataFrame([
        {"sample_id": "EC001", "expected_genes": "blaCTX-M-15"},
        {"sample_id": "EC002", "expected_genes": "none"},
        {"sample_id": "EC003", "expected_genes": "blaCMY-2"},
    ])

    fastp_table = pd.DataFrame([
        {"sample_id": "EC001", "gc_content_percent": 50.1, "coverage_status": "PASS"},
        {"sample_id": "EC002", "gc_content_percent": 49.8, "coverage_status": "PASS"},
        {"sample_id": "EC003", "gc_content_percent": 50.4, "coverage_status": "WARNING"},
    ])
    quast_table = pd.DataFrame([
        {"sample_id": "EC001", "gc_content_percent": 50.65, "assembly_status": "PASS"},
        {"sample_id": "EC002", "gc_content_percent": 50.20, "assembly_status": "PASS"},
        {"sample_id": "EC003", "gc_content_percent": 50.40, "assembly_status": "PASS"},
    ])
    checkm_table = pd.DataFrame([
        {"sample_id": "EC001", "completeness_status": "PASS"},
        {"sample_id": "EC002", "completeness_status": "FAIL"},  # contaminacion alta
        {"sample_id": "EC003", "completeness_status": "PASS"},
    ])
    taxonomy_table = pd.DataFrame([
        {"sample_id": "EC001", "taxonomy_status": "PASS"},
        {"sample_id": "EC002", "taxonomy_status": "PASS"},
        {"sample_id": "EC003", "taxonomy_status": "PASS"},
    ])
    amr_table = pd.DataFrame([
        {"sample_id": "EC001", "gene_symbol": "blaCTX-M-15", "antimicrobial_class": "BETA-LACTAM", "meets_identity_coverage_threshold": True},
        {"sample_id": "EC003", "gene_symbol": "blaCMY-2", "antimicrobial_class": "BETA-LACTAM", "meets_identity_coverage_threshold": True},
    ])
    reference_comparison_table = pd.DataFrame([
        {"sample_id": "EC001", "detected_gene": "blaCTX-M-15", "match_type": "exact", "pipeline_status": "positive", "confusion_category": "TP"},
        {"sample_id": "EC002", "detected_gene": "none", "match_type": "none", "pipeline_status": "negative", "confusion_category": "TN"},
        {"sample_id": "EC003", "detected_gene": "blaCMY-2", "match_type": "exact", "pipeline_status": "positive", "confusion_category": "TP"},
    ])
    performance_table = pd.DataFrame([
        {"sample_id": "EC001", "total_elapsed_seconds": 661.1, "peak_max_ram_gb": 7.81},
        {"sample_id": "EC002", "total_elapsed_seconds": 615.3, "peak_max_ram_gb": 7.20},
        {"sample_id": "EC003", "total_elapsed_seconds": 19.0, "peak_max_ram_gb": 0.41},
    ])

    master_table = build_master_table(
        samples_table, fastp_table, quast_table, checkm_table, taxonomy_table,
        amr_table, reference_comparison_table, performance_table,
    ).set_index("sample_id")

    # EC001: todo PASS -> PASS.
    assert master_table.loc["EC001", "final_status"] == "PASS"
    assert master_table.loc["EC001", "detected_gene_count"] == 1
    # EC002: CheckM FAIL -> EXCLUDED, pero sigue apareciendo en la tabla.
    assert master_table.loc["EC002", "final_status"] == "EXCLUDED"
    # EC003: fastp WARNING, resto PASS -> WARNING.
    assert master_table.loc["EC003", "final_status"] == "WARNING"

    # Las columnas "gc_content_percent" de fastp y QUAST no deben colisionar.
    assert master_table.loc["EC001", "reads_gc_content_percent"] == 50.1
    assert master_table.loc["EC001", "assembly_gc_content_percent"] == 50.65

    # El tiempo/RAM de desempeno llega correctamente a la tabla maestra.
    assert master_table.loc["EC001", "total_elapsed_seconds"] == 661.1
    assert master_table.loc["EC001", "peak_max_ram_gb"] == 7.81


def test_master_table_tolerates_missing_module_tables():
    # Corrida parcial durante el desarrollo: si un modulo todavia no corrio
    # para ninguna muestra, la tabla maestra se arma igual con lo disponible.
    samples_table = pd.DataFrame([{"sample_id": "EC001", "expected_genes": "none"}])

    master_table = build_master_table(
        samples_table,
        fastp_table=None, quast_table=None, checkm_table=None, taxonomy_table=None,
        amr_table=None, reference_comparison_table=None, performance_by_sample_table=None,
    )

    assert list(master_table["sample_id"]) == ["EC001"]
    assert master_table.iloc[0]["final_status"] == "PENDING"
