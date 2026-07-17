"""Concordancia analitica entre dos motores independientes de deteccion de AMR
(AMRFinderPlus y ABricate/CARD+ResFinder).

Este es un pilar distinto de la comparacion contra el estandar de referencia
fenotipico (compare_to_reference.py, seccion 14 del diseno original) y de la
reproducibilidad entre corridas repetidas (assess_reproducibility.py): aqui se
mide si DOS HERRAMIENTAS DE DETECCION INDEPENDIENTES, corridas sobre el MISMO
ensamblaje, coinciden en que genes de resistencia encuentran. Una concordancia
baja entre motores es una senal de alerta metodologica -- distinta de una
discrepancia contra el fenotipo, que podria deberse a expresion genica, no a
un error de deteccion.

La comparacion se hace a nivel de FAMILIA de gen (columna gene_family, ya
derivada por parse_amrfinder.py y parse_abricate.py), nunca de alelo exacto:
AMRFinderPlus (catalogo de referencia de NCBI) y CARD/ResFinder (via ABricate)
no siempre nombran el mismo alelo de forma identica, y comparar alelo por
alelo subestimaria artificialmente la concordancia real entre herramientas.

Produce dos salidas:

  results/tables/engine_concordance.tsv
      Una fila por muestra: concordancia exacta (Jaccard = 1.0 o no) y
      similitud de Jaccard entre los conjuntos de familias detectadas por
      cada motor. Nunca se calcula ningun estadistico formal aqui (eso es
      exclusivo de R, ver compare_engines_statistics.R).

  results/statistics/engine_concordance_input.csv
      Tabla larga (una fila por combinacion muestra+familia detectada por AL
      MENOS un motor), con el resultado de cada motor como factor
      detected/not_detected -- lista para que R calcule el indice kappa
      entre motores.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

ENGINE_CONCORDANCE_COLUMNS = [
    "sample_id",
    "amrfinder_gene_families",
    "abricate_gene_families",
    "exact_gene_family_concordance",
    "jaccard_similarity",
]

ENGINE_AGREEMENT_COLUMNS = ["sample_id", "gene_family", "amrfinder_result", "abricate_result"]


def exact_gene_concordance(families_engine_a: set[str], families_engine_b: set[str]) -> float:
    """1.0 si dos motores detectaron exactamente el mismo conjunto de
    familias de gen, 0.0 si no (misma logica que assess_reproducibility.py,
    aplicada aqui entre motores en vez de entre corridas repetidas)."""
    return 1.0 if families_engine_a == families_engine_b else 0.0


def jaccard_similarity(families_engine_a: set[str], families_engine_b: set[str]) -> float:
    """Fraccion de familias compartidas sobre el total de familias distintas
    detectadas por cualquiera de los dos motores."""
    union = families_engine_a | families_engine_b
    if not union:
        return 1.0
    return len(families_engine_a & families_engine_b) / len(union)


def build_gene_family_sets_by_sample(normalized_long_table: pd.DataFrame) -> dict[str, set[str]]:
    """Conjunto de familias de gen detectadas con confianza, por muestra, a
    partir de una tabla larga ya normalizada (parse_amrfinder.py o
    parse_abricate.py, ambas comparten el mismo esquema de columnas)."""
    confident_detections = normalized_long_table.loc[normalized_long_table["meets_identity_coverage_threshold"]]
    return confident_detections.groupby("sample_id")["gene_family"].apply(set).to_dict()


def build_engine_concordance_report(
    sample_ids: list[str],
    amrfinder_families_by_sample: dict[str, set[str]],
    abricate_families_by_sample: dict[str, set[str]],
) -> pd.DataFrame:
    """Arma una fila por muestra comparando los conjuntos de familias de gen
    detectados por cada motor. Toda muestra documentada aparece, incluso si
    ningun motor detecto nada en ella."""
    concordance_rows = []
    for sample_id in sample_ids:
        families_amrfinder = amrfinder_families_by_sample.get(sample_id, set())
        families_abricate = abricate_families_by_sample.get(sample_id, set())
        concordance_rows.append({
            "sample_id": sample_id,
            "amrfinder_gene_families": ", ".join(sorted(families_amrfinder)) if families_amrfinder else "none",
            "abricate_gene_families": ", ".join(sorted(families_abricate)) if families_abricate else "none",
            "exact_gene_family_concordance": exact_gene_concordance(families_amrfinder, families_abricate),
            "jaccard_similarity": round(jaccard_similarity(families_amrfinder, families_abricate), 4),
        })
    return pd.DataFrame(concordance_rows, columns=ENGINE_CONCORDANCE_COLUMNS)


def build_engine_agreement_long_table(
    sample_ids: list[str],
    amrfinder_families_by_sample: dict[str, set[str]],
    abricate_families_by_sample: dict[str, set[str]],
) -> pd.DataFrame:
    """Arma la tabla larga (muestra x familia de gen) que R usara para
    calcular el indice kappa entre motores. Solo incluye familias detectadas
    por AL MENOS un motor en esa muestra -- una familia que ningun motor
    detecto en ninguna muestra no aporta informacion a la concordancia."""
    agreement_rows = []
    for sample_id in sample_ids:
        families_amrfinder = amrfinder_families_by_sample.get(sample_id, set())
        families_abricate = abricate_families_by_sample.get(sample_id, set())
        for gene_family in sorted(families_amrfinder | families_abricate):
            agreement_rows.append({
                "sample_id": sample_id,
                "gene_family": gene_family,
                "amrfinder_result": "detected" if gene_family in families_amrfinder else "not_detected",
                "abricate_result": "detected" if gene_family in families_abricate else "not_detected",
            })
    return pd.DataFrame(agreement_rows, columns=ENGINE_AGREEMENT_COLUMNS)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comparar los genes detectados por AMRFinderPlus y ABricate (concordancia entre motores)."
    )
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument(
        "--amrfinder-table", type=Path, default=Path("results/tables/amr_summary.tsv"),
        help="Listado largo normalizado de AMRFinderPlus (parse_amrfinder.py combine)",
    )
    parser.add_argument(
        "--abricate-table", type=Path, default=Path("results/tables/abricate_summary.tsv"),
        help="Listado largo normalizado de ABricate (parse_abricate.py combine)",
    )
    parser.add_argument(
        "--concordance-output", type=Path, default=Path("results/tables/engine_concordance.tsv"),
    )
    parser.add_argument(
        "--agreement-output", type=Path, default=Path("results/statistics/engine_concordance_input.csv"),
    )
    args = parser.parse_args()

    samples_table = pd.read_csv(args.samples, sep="\t", dtype=str)
    sample_ids = samples_table["sample_id"].tolist()

    amrfinder_table = pd.read_csv(args.amrfinder_table, sep="\t")
    abricate_table = pd.read_csv(args.abricate_table, sep="\t")

    amrfinder_families_by_sample = build_gene_family_sets_by_sample(amrfinder_table)
    abricate_families_by_sample = build_gene_family_sets_by_sample(abricate_table)

    concordance_report = build_engine_concordance_report(
        sample_ids, amrfinder_families_by_sample, abricate_families_by_sample
    )
    args.concordance_output.parent.mkdir(parents=True, exist_ok=True)
    concordance_report.to_csv(args.concordance_output, sep="\t", index=False)

    agreement_table = build_engine_agreement_long_table(
        sample_ids, amrfinder_families_by_sample, abricate_families_by_sample
    )
    args.agreement_output.parent.mkdir(parents=True, exist_ok=True)
    agreement_table.to_csv(args.agreement_output, index=False)

    fully_concordant = (concordance_report["exact_gene_family_concordance"] == 1.0).sum()
    print(
        f"{len(concordance_report)} muestra(s) comparadas entre motores, "
        f"{fully_concordant} totalmente concordantes, escritas en {args.concordance_output}"
    )
    print(f"{len(agreement_table)} fila(s) de acuerdo muestra+familia escritas en {args.agreement_output}")


if __name__ == "__main__":
    main()
