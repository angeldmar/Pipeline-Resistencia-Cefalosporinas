"""Comparacion de los genes detectados por el pipeline contra el estandar de
referencia (el gen esperado que viene documentado en config/samples.tsv).

Este script arma la matriz de datos (TP/TN/FP/FN) que R usara mas adelante
para calcular sensibilidad, especificidad, kappa, etc. Python solo prepara y
clasifica los datos; NINGUN calculo estadistico se hace aqui (eso es
exclusivo de R, ver seccion de estadistica del diseno del pipeline).

Alcance de esta comparacion: es puramente a nivel de gen (lo que detecto
AMRFinderPlus vs. lo que se esperaba segun expected_genes). Todavia NO tiene
en cuenta si la muestra fallo algun control de calidad anterior (cobertura,
ensamblaje, taxonomia, completitud) -- esa integracion ocurre en la tabla
maestra (merge_results.py), que es quien tiene visibilidad de todos los
modulos a la vez. Aqui "Indeterminado" se usa unicamente cuando el propio
dato de referencia (expected_genes) esta vacio o marcado como NA, es decir,
cuando no hay estandar de referencia contra el cual comparar.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re

import pandas as pd

REFERENCE_COMPARISON_COLUMNS = [
    "sample_id",
    "expected_gene",
    "detected_gene",
    "match_type",
    "reference_status",
    "pipeline_status",
    "confusion_category",
]

# Valores de expected_genes que significan "no se espera ningun gen de
# resistencia" (referencia negativa), distintos de un dato de referencia
# faltante/desconocido (ver NO_REFERENCE_VALUES abajo).
NEGATIVE_REFERENCE_VALUES = {"none"}
# Valores que significan "no hay estandar de referencia documentado para
# esta muestra" -> no se puede clasificar como TP/TN/FP/FN, solo "Indeterminado".
NO_REFERENCE_VALUES = {"", "na", "nan"}

# Misma heuristica de familia de gen que parse_amrfinder.py (se duplica aqui
# a proposito: cada script de este pipeline es autocontenido y se puede
# invocar de forma independiente, sin depender de importar otros scripts).
ALLELE_SUFFIX_PATTERN = re.compile(r"[-_]\d+(\.\d+)?[A-Za-z]?$")


def derive_gene_family(gene_symbol: str) -> str:
    return ALLELE_SUFFIX_PATTERN.sub("", gene_symbol)


def determine_reference_status(expected_gene: str) -> tuple[bool | None, str]:
    """Interpreta el valor de expected_genes de una muestra.

    Devuelve (True, "positive") si se espera un gen de resistencia,
    (False, "negative") si explicitamente no se espera ninguno ("none"), o
    (None, "indeterminate") si no hay dato de referencia documentado.
    """
    normalized_value = str(expected_gene).strip().lower()
    if normalized_value in NO_REFERENCE_VALUES:
        return None, "indeterminate"
    if normalized_value in NEGATIVE_REFERENCE_VALUES:
        return False, "negative"
    return True, "positive"


def confusion_category(reference_positive: bool | None, pipeline_positive: bool) -> str:
    """Clasifica la comparacion en TP/TN/FP/FN, o "Indeterminado" cuando no
    hay estandar de referencia contra el cual comparar (reference_positive
    es None). Para los cuatro casos con referencia conocida, esta funcion es
    identica a la del diseno del pipeline."""
    if reference_positive is None:
        return "Indeterminado"
    if reference_positive and pipeline_positive:
        return "TP"
    if not reference_positive and not pipeline_positive:
        return "TN"
    if not reference_positive and pipeline_positive:
        return "FP"
    return "FN"


def build_comparison_row(sample_id: str, expected_gene: str, sample_amr_rows: pd.DataFrame) -> dict:
    """Compara el gen esperado de una muestra contra sus genes detectados
    (ya filtrados a los que cumplen el umbral de identidad/cobertura)."""
    reference_positive, reference_status = determine_reference_status(expected_gene)

    detected_gene_symbols = set(sample_amr_rows["gene_symbol"])
    detected_families_to_symbol = {derive_gene_family(symbol): symbol for symbol in detected_gene_symbols}

    # reference_positive es un bool O None (indeterminado): se comparan las
    # tres ramas de forma explicita en vez de usar "if reference_positive",
    # porque None tambien evalua como falso en Python y se confundiria
    # silenciosamente con el caso "referencia negativa".
    if reference_positive is None:
        # Sin estandar de referencia documentado: se informa igual lo que el
        # pipeline detecto (si algo), pero sin calificarlo de "esperado" ni
        # "inesperado" -- no hay nada contra que comparar. La categoria final
        # queda en "Indeterminado" sin importar pipeline_positive.
        detected_beta_lactam_genes = sorted(
            sample_amr_rows.loc[
                sample_amr_rows["antimicrobial_class"].astype(str).str.upper() == "BETA-LACTAM",
                "gene_symbol",
            ].unique()
        )
        detected_gene = ", ".join(detected_beta_lactam_genes) if detected_beta_lactam_genes else "none"
        match_type = "not_applicable"
        pipeline_positive = bool(detected_beta_lactam_genes)
    elif reference_positive:
        if expected_gene in detected_gene_symbols:
            detected_gene, match_type, pipeline_positive = expected_gene, "exact", True
        elif derive_gene_family(expected_gene) in detected_families_to_symbol:
            # Mismo gen, alelo distinto al documentado como referencia (ej.
            # se esperaba blaCTX-M-15 y se detecto blaCTX-M-27): se cuenta
            # como deteccion positiva, mismo criterio del gen pero se deja
            # constancia de que el alelo exacto no coincidio.
            detected_gene = detected_families_to_symbol[derive_gene_family(expected_gene)]
            match_type, pipeline_positive = "family", True
        else:
            detected_gene, match_type, pipeline_positive = "none", "none", False
    else:
        # No se espera ningun gen: cualquier beta-lactamasa detectada con
        # confianza suficiente es una posible discrepancia con la referencia.
        unexpected_beta_lactam_genes = sorted(
            sample_amr_rows.loc[
                sample_amr_rows["antimicrobial_class"].astype(str).str.upper() == "BETA-LACTAM",
                "gene_symbol",
            ].unique()
        )
        if unexpected_beta_lactam_genes:
            detected_gene = ", ".join(unexpected_beta_lactam_genes)
            match_type, pipeline_positive = "unexpected", True
        else:
            detected_gene, match_type, pipeline_positive = "none", "none", False

    pipeline_status = "positive" if pipeline_positive else "negative"
    category = confusion_category(reference_positive, pipeline_positive)

    return {
        "sample_id": sample_id,
        "expected_gene": expected_gene,
        "detected_gene": detected_gene,
        "match_type": match_type,
        "reference_status": reference_status,
        "pipeline_status": pipeline_status,
        "confusion_category": category,
    }


def build_comparison_table(samples_table: pd.DataFrame, amr_table: pd.DataFrame) -> pd.DataFrame:
    """Construye la tabla de comparacion completa, una fila por muestra."""
    confident_amr_table = amr_table.loc[amr_table["meets_identity_coverage_threshold"]]

    comparison_rows = []
    for _, sample_row in samples_table.iterrows():
        sample_id = sample_row["sample_id"]
        sample_amr_rows = confident_amr_table.loc[confident_amr_table["sample_id"] == sample_id]
        comparison_rows.append(
            build_comparison_row(sample_id, sample_row["expected_genes"], sample_amr_rows)
        )

    return pd.DataFrame(comparison_rows, columns=REFERENCE_COMPARISON_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comparar los genes detectados por el pipeline contra el estandar de referencia."
    )
    parser.add_argument(
        "--samples", type=Path, default=Path("config/samples.tsv"),
        help="Tabla de muestras validada, con la columna expected_genes",
    )
    parser.add_argument(
        "--amr-table", type=Path, default=Path("results/tables/amr_summary.tsv"),
        help="Listado largo de AMR ya normalizado (results/tables/amr_summary.tsv)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/tables/reference_comparison.tsv"),
        help="Ruta de salida de la matriz de comparacion",
    )
    args = parser.parse_args()

    # dtype=str + fillna("NA"): sin esto, pandas interpreta el texto literal
    # "NA" de la columna expected_genes como valor nulo (NaN) en vez de como
    # el string "NA", igual que ya se maneja en validate_samples.py.
    samples_table = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("NA")
    amr_table = pd.read_csv(args.amr_table, sep="\t")

    comparison_table = build_comparison_table(samples_table, amr_table)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    comparison_table.to_csv(args.output, sep="\t", index=False)

    category_counts = comparison_table["confusion_category"].value_counts().to_dict()
    print(f"{len(comparison_table)} muestra(s) comparadas, escritas en {args.output}")
    print(f"Distribucion: {category_counts}")


if __name__ == "__main__":
    main()
