"""Integracion de resultados: arma la tabla maestra del pipeline.

Junta, en una sola tabla con UNA FILA POR MUESTRA, los resultados de todos
los modulos anteriores: metadatos y fenotipo (samples.tsv), calidad de
lecturas y cobertura (fastp), metricas de ensamblaje (QUAST), completitud y
contaminacion (CheckM), identificacion taxonomica (Kraken2), un resumen de
los genes de resistencia detectados (AMRFinderPlus), la comparacion contra
el estandar de referencia, y el tiempo/RAM totales por muestra.

Siguiendo la recomendacion de la seccion de integracion de resultados del
diseno del pipeline, esta tabla maestra NO intenta meter el detalle de cada
gen detectado en columnas (eso generaria una tabla ancha y dificil de leer,
distinta por muestra segun cuantos genes tenga). El detalle gen por gen ya
vive en la tabla larga (results/tables/amr_summary.tsv /
amr_classified.tsv); aqui solo se agrega un resumen corto (cuantos genes,
cuales) para tener una vista rapida por muestra. Lo mismo aplica al
desempeno computacional: aqui solo entran el tiempo TOTAL y la RAM maxima
por muestra; el detalle por modulo vive en performance_summary.tsv /
performance_by_module.tsv (ver combine_performance.py). La concordancia
entre motores de AMR (AMRFinderPlus vs. ABricate, ver compare_amr_engines.py)
tambien se resume en dos columnas cortas; el detalle por familia de gen vive
en engine_concordance.tsv.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd

# Nombres de columna que colisionan entre modulos (ej. tanto fastp como QUAST
# reportan "gc_content_percent", pero de cosas distintas: lecturas crudas vs.
# ensamblaje) y necesitan renombrarse antes de combinarlos en la tabla maestra.
FASTP_COLUMN_RENAMES = {"gc_content_percent": "reads_gc_content_percent"}
QUAST_COLUMN_RENAMES = {"gc_content_percent": "assembly_gc_content_percent"}

# Los estados de modulo (PASS/WARNING/FAIL) que se consideran para decidir el
# estado final combinado de una muestra.
QC_GATE_STATUS_COLUMNS = ["coverage_status", "assembly_status", "completeness_status", "taxonomy_status"]


def load_optional_table(table_path: Path, module_name: str) -> pd.DataFrame | None:
    """Carga la tabla de resumen de un modulo si existe. Si todavia no se ha
    generado (por ejemplo, porque ese modulo no se ha corrido para ninguna
    muestra), se avisa por stdout y se continua sin el en vez de fallar: la
    tabla maestra se arma con lo que si este disponible."""
    if not table_path.is_file():
        print(f"AVISO: no se encontro la tabla de {module_name} en {table_path}; se omite en la tabla maestra.")
        return None
    return pd.read_csv(table_path, sep="\t")


def build_amr_sample_summary(amr_table: pd.DataFrame) -> pd.DataFrame:
    """Reduce el listado largo de AMR (una fila por gen) a un resumen corto
    por muestra: cuantos genes confiables se detectaron en total, cuantos son
    de la clase BETA-LACTAM, y su listado (para una vista rapida; el detalle
    completo -- identidad, cobertura, coordenadas -- sigue viviendo en
    amr_summary.tsv / amr_classified.tsv)."""
    confident_detections = amr_table.loc[amr_table["meets_identity_coverage_threshold"]]

    summary_columns = ["sample_id", "detected_gene_count", "detected_beta_lactam_gene_count", "detected_genes"]
    if confident_detections.empty:
        return pd.DataFrame(columns=summary_columns)

    gene_counts_and_lists = confident_detections.groupby("sample_id").agg(
        detected_gene_count=("gene_symbol", "count"),
        detected_genes=("gene_symbol", lambda gene_symbols: ", ".join(sorted(gene_symbols))),
    ).reset_index()

    beta_lactam_detections = confident_detections.loc[
        confident_detections["antimicrobial_class"].astype(str).str.upper() == "BETA-LACTAM"
    ]
    beta_lactam_counts = (
        beta_lactam_detections.groupby("sample_id").size().rename("detected_beta_lactam_gene_count").reset_index()
    )

    combined_summary = gene_counts_and_lists.merge(beta_lactam_counts, on="sample_id", how="left")
    return combined_summary[summary_columns]


def determine_final_status(qc_gate_statuses: list[str]) -> str:
    """Combina los estados de los modulos de control de calidad (cobertura,
    ensamblaje, completitud/contaminacion, taxonomia) en un veredicto final.

    Una muestra nunca se elimina silenciosamente: si algun modulo dio FAIL,
    el estado final es EXCLUDED (candidata a excluirse del analisis
    principal, pero permanece en la tabla, visible, con el motivo
    identificable en sus columnas de origen). WARNING si ningun modulo fallo
    pero al menos uno advirtio. PASS si todos los modulos disponibles
    pasaron. PENDING si todavia no hay ningun modulo de QC corrido para esa
    muestra (por ejemplo, en una corrida parcial durante el desarrollo).
    """
    known_statuses = [status for status in qc_gate_statuses if pd.notna(status)]
    if not known_statuses:
        return "PENDING"
    if "FAIL" in known_statuses:
        return "EXCLUDED"
    if "WARNING" in known_statuses:
        return "WARNING"
    return "PASS"


def build_master_table(
    samples_table: pd.DataFrame,
    fastp_table: pd.DataFrame | None,
    quast_table: pd.DataFrame | None,
    checkm_table: pd.DataFrame | None,
    taxonomy_table: pd.DataFrame | None,
    amr_table: pd.DataFrame | None,
    reference_comparison_table: pd.DataFrame | None,
    performance_by_sample_table: pd.DataFrame | None = None,
    engine_concordance_table: pd.DataFrame | None = None,
    mlst_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Arma la tabla maestra: parte de samples.tsv (para que TODA muestra
    documentada aparezca, incluso si a algun modulo todavia le falta
    correr) y le va uniendo (left join, por sample_id) el resumen de cada
    modulo disponible."""
    master_table = samples_table.copy()

    if fastp_table is not None:
        master_table = master_table.merge(
            fastp_table.rename(columns=FASTP_COLUMN_RENAMES), on="sample_id", how="left"
        )

    if quast_table is not None:
        master_table = master_table.merge(
            quast_table.rename(columns=QUAST_COLUMN_RENAMES), on="sample_id", how="left"
        )

    if checkm_table is not None:
        master_table = master_table.merge(checkm_table, on="sample_id", how="left")

    if taxonomy_table is not None:
        master_table = master_table.merge(taxonomy_table, on="sample_id", how="left")

    if amr_table is not None:
        amr_sample_summary = build_amr_sample_summary(amr_table)
        master_table = master_table.merge(amr_sample_summary, on="sample_id", how="left")
        master_table["detected_gene_count"] = master_table["detected_gene_count"].fillna(0).astype(int)
        master_table["detected_beta_lactam_gene_count"] = (
            master_table["detected_beta_lactam_gene_count"].fillna(0).astype(int)
        )
        master_table["detected_genes"] = master_table["detected_genes"].fillna("none")

    if reference_comparison_table is not None:
        reference_columns = reference_comparison_table[
            ["sample_id", "detected_gene", "match_type", "pipeline_status", "confusion_category"]
        ].rename(columns={"detected_gene": "reference_detected_gene"})
        master_table = master_table.merge(reference_columns, on="sample_id", how="left")

    if performance_by_sample_table is not None:
        # Solo el tiempo TOTAL y la RAM MAXIMA por muestra entran a la tabla
        # maestra; el desglose por modulo vive en performance_summary.tsv /
        # performance_by_module.tsv (evita repetir aqui el mismo problema de
        # "columnas dificiles de analizar" que ya se evito con los genes de AMR).
        performance_columns = performance_by_sample_table[
            ["sample_id", "total_elapsed_seconds", "peak_max_ram_gb"]
        ]
        master_table = master_table.merge(performance_columns, on="sample_id", how="left")

    if engine_concordance_table is not None:
        # Incluye tambien los conjuntos de familias detectados por cada
        # motor (texto corto, util para interpretar una discordancia sin
        # tener que abrir engine_concordance.tsv aparte).
        engine_columns = engine_concordance_table[
            ["sample_id", "amrfinder_gene_families", "abricate_gene_families",
             "exact_gene_family_concordance", "jaccard_similarity"]
        ].rename(columns={
            "exact_gene_family_concordance": "engine_exact_concordance",
            "jaccard_similarity": "engine_jaccard_similarity",
        })
        master_table = master_table.merge(engine_columns, on="sample_id", how="left")

    if mlst_table is not None:
        # Solo el ST y su confiabilidad; el perfil de alelos completo
        # (allele_profile) queda en mlst_summary.tsv para no ensanchar la
        # tabla maestra con un texto largo por muestra.
        mlst_columns = mlst_table[["sample_id", "scheme", "sequence_type", "sequence_type_status"]]
        master_table = master_table.merge(mlst_columns, on="sample_id", how="left")

    available_qc_gate_columns = [
        column for column in QC_GATE_STATUS_COLUMNS if column in master_table.columns
    ]
    master_table["final_status"] = master_table[available_qc_gate_columns].apply(
        lambda row: determine_final_status(list(row)), axis=1
    )

    return master_table


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Combinar los resultados de todos los modulos del pipeline en una tabla maestra."
    )
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument("--fastp-summary", type=Path, default=Path("results/tables/fastp_summary.tsv"))
    parser.add_argument("--quast-summary", type=Path, default=Path("results/tables/quast_summary.tsv"))
    parser.add_argument("--checkm-summary", type=Path, default=Path("results/tables/checkm_summary.tsv"))
    parser.add_argument("--taxonomy-summary", type=Path, default=Path("results/tables/taxonomy_summary.tsv"))
    parser.add_argument("--amr-summary", type=Path, default=Path("results/tables/amr_summary.tsv"))
    parser.add_argument(
        "--reference-comparison", type=Path, default=Path("results/tables/reference_comparison.tsv")
    )
    parser.add_argument(
        "--performance-by-sample", type=Path, default=Path("results/tables/performance_by_sample.tsv"),
        help="Resumen de tiempo total / RAM maxima por muestra (ver combine_performance.py)",
    )
    parser.add_argument(
        "--engine-concordance", type=Path, default=Path("results/tables/engine_concordance.tsv"),
        help="Concordancia AMRFinderPlus vs. ABricate por muestra (ver compare_amr_engines.py)",
    )
    parser.add_argument(
        "--mlst-summary", type=Path, default=Path("results/tables/mlst_summary.tsv"),
        help="Tipificacion de secuencia multilocus por muestra (ver parse_mlst.py)",
    )
    parser.add_argument("--output", type=Path, default=Path("results/tables/master_results.tsv"))
    args = parser.parse_args()

    samples_table = pd.read_csv(args.samples, sep="\t", dtype=str).fillna("NA")

    master_table = build_master_table(
        samples_table,
        load_optional_table(args.fastp_summary, "calidad/cobertura (fastp)"),
        load_optional_table(args.quast_summary, "ensamblaje (QUAST)"),
        load_optional_table(args.checkm_summary, "completitud/contaminacion (CheckM)"),
        load_optional_table(args.taxonomy_summary, "taxonomia (Kraken2)"),
        load_optional_table(args.amr_summary, "deteccion de AMR (AMRFinderPlus)"),
        load_optional_table(args.reference_comparison, "comparacion con el estandar de referencia"),
        load_optional_table(args.performance_by_sample, "desempeño computacional"),
        load_optional_table(args.engine_concordance, "concordancia entre motores de AMR"),
        load_optional_table(args.mlst_summary, "tipificacion MLST"),
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    master_table.to_csv(args.output, sep="\t", index=False)

    final_status_counts = master_table["final_status"].value_counts().to_dict()
    print(f"Tabla maestra de {len(master_table)} muestra(s) escrita en {args.output}")
    print(f"Estado final: {final_status_counts}")


if __name__ == "__main__":
    main()
