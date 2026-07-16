"""Prueba de integracion: normalizacion -> clasificacion -> comparacion de AMR.

Alcance: no hay forma de correr AMRFinderPlus real en este entorno de
desarrollo (ver seccion "ambientes Conda", QUAST/CheckM/Prokka ni siquiera
tienen build para esta arquitectura). Esta prueba verifica en cambio el
contrato REAL entre nuestros propios scripts -- que es donde vive el riesgo
real de integracion en este pipeline (los formatos que un script produce y
el siguiente consume) -- encadenando parse_amrfinder.py ->
classify_cephalosporin_genes.py -> compare_to_reference.py con datos que
imitan el formato real de salida de AMRFinderPlus.
"""

import pandas as pd

from classify_cephalosporin_genes import classify_amr_table, load_resistance_targets
from compare_to_reference import build_comparison_table
from parse_amrfinder import normalize_amrfinder_table


def build_raw_amrfinder_table() -> pd.DataFrame:
    """Imita el formato de columnas real de AMRFinderPlus (modo --nucleotide)."""
    return pd.DataFrame([
        {
            "Gene symbol": "blaCTX-M-15", "Sequence name": "CTX-M-15 family class A ESBL",
            "Class": "BETA-LACTAM", "Subclass": "CEPHALOSPORIN", "Method": "ALLELEX",
            "% Identity to reference sequence": 100.0, "% Coverage of reference sequence": 100.0,
            "Contig id": "contig_1", "Start": 1, "Stop": 900,
        },
        {
            "Gene symbol": "blaTEM-1", "Sequence name": "TEM family class A beta-lactamase",
            "Class": "BETA-LACTAM", "Subclass": "BETA-LACTAM", "Method": "ALLELEX",
            "% Identity to reference sequence": 99.5, "% Coverage of reference sequence": 100.0,
            "Contig id": "contig_2", "Start": 500, "Stop": 1360,
        },
    ])


def test_normalized_amr_table_flows_correctly_into_classification_and_comparison(repo_root):
    raw_table = build_raw_amrfinder_table()

    normalized_table = normalize_amrfinder_table(
        "EC001", raw_table, minimum_identity=90, minimum_gene_coverage=80
    )

    # El formato normalizado debe tener las columnas que classify_amr_table
    # y build_comparison_table esperan, sin ningun paso de traduccion manual.
    resistance_targets = load_resistance_targets(repo_root / "config" / "resistance_targets.yaml")
    classified_table = classify_amr_table(normalized_table, resistance_targets)

    ctx_m_row = classified_table.loc[classified_table["gene_symbol"] == "blaCTX-M-15"].iloc[0]
    tem_row = classified_table.loc[classified_table["gene_symbol"] == "blaTEM-1"].iloc[0]
    assert ctx_m_row["beta_lactamase_category"] == "ESBL"
    assert tem_row["beta_lactamase_category"] == "Other"  # no se asume BLEE solo por el prefijo

    # La misma tabla normalizada (no la clasificada) es la que consume
    # compare_to_reference.py -- se verifica que tambien encaje ahi.
    samples_table = pd.DataFrame([{"sample_id": "EC001", "expected_genes": "blaCTX-M-15"}])
    comparison_table = build_comparison_table(samples_table, normalized_table)

    assert comparison_table.iloc[0]["confusion_category"] == "TP"


def test_low_confidence_detection_excluded_from_reference_comparison_but_visible_in_normalized_table():
    # Un unico gen, con cobertura por debajo del umbral (50% < 80%).
    raw_table = build_raw_amrfinder_table().iloc[[1]].copy()
    raw_table["% Coverage of reference sequence"] = 50.0

    normalized_table = normalize_amrfinder_table(
        "EC001", raw_table, minimum_identity=90, minimum_gene_coverage=80
    )

    # El gen de baja confianza sigue en la tabla normalizada (no se descarta)...
    assert len(normalized_table) == 1
    assert not normalized_table.iloc[0]["meets_identity_coverage_threshold"]

    # ...pero compare_to_reference.py, que solo mira detecciones confiables,
    # no lo cuenta como una deteccion valida.
    samples_table = pd.DataFrame([{"sample_id": "EC001", "expected_genes": "none"}])
    comparison_table = build_comparison_table(samples_table, normalized_table)
    assert comparison_table.iloc[0]["confusion_category"] == "TN"  # no FP, porque el gen no cuenta
