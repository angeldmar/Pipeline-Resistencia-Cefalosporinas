"""Pruebas de reproducibilidad: compara corridas repetidas de la misma muestra.

Cada muestra de la evaluacion de reproducibilidad se procesa al menos tres
veces bajo las mismas condiciones, con sample_id del estilo "EC001_run1",
"EC001_run2", "EC001_run3" (mismo aislamiento, tres corridas independientes
del pipeline). Este script agrupa esas corridas por su muestra base y
compara, entre cada par de corridas:

  - Genes detectados y sus alelos (concordancia exacta + Jaccard).
  - El archivo de ensamblaje final (hash SHA-256: identico o no).
  - El estado de clasificacion final (final_status de la tabla maestra).

IMPORTANTE sobre el coeficiente de variacion: el diseno del pipeline reserva
el calculo de CV exclusivamente para R ("R se reservara unicamente para...
Coeficiente de variacion"). Por eso este script NO calcula CV, ni siquiera
como medida secundaria -- solo usa comparaciones categoricas/exactas
(concordancia exacta, Jaccard, igualdad de hash). Los datos numericos
crudos por corrida (cobertura, tiempo, RAM) quedan disponibles en
master_results.tsv con su columna "run", listos para que R calcule el CV
formalmente en el paso de estadistica.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import itertools
import re

import pandas as pd

REPLICATE_RUN_ID_PATTERN = re.compile(r"^(?P<base_sample_id>.+)_run(?P<run_number>\d+)$")


def parse_replicate_run_id(sample_id: str) -> tuple[str, int] | None:
    """Separa un sample_id de reproducibilidad (ej. "EC001_run2") en su
    muestra base ("EC001") y numero de corrida (2). Devuelve None si el
    sample_id no sigue esa convencion (ej. una muestra normal, sin repetir)."""
    match = REPLICATE_RUN_ID_PATTERN.match(sample_id)
    if not match:
        return None
    return match.group("base_sample_id"), int(match.group("run_number"))


def exact_gene_concordance(genes_run_a: set[str], genes_run_b: set[str]) -> float:
    """1.0 si dos corridas detectaron exactamente el mismo conjunto de genes
    (mismos simbolos, que en AMRFinderPlus ya incluyen el alelo), 0.0 si no."""
    return 1.0 if genes_run_a == genes_run_b else 0.0


def jaccard_similarity(set_a: set[str], set_b: set[str]) -> float:
    """Fraccion de genes compartidos sobre el total de genes distintos entre
    las dos corridas. Mas informativo que la concordancia exacta cuando dos
    corridas casi coinciden pero difieren en un solo gen."""
    union = set_a | set_b
    if not union:
        return 1.0
    return len(set_a & set_b) / len(union)


def compute_sha256(file_path: Path) -> str:
    """Calcula el hash SHA-256 de un archivo, leyendolo en bloques (misma
    logica que download_data.py; se duplica aqui a proposito, ya que cada
    script de este pipeline es autocontenido)."""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as file_handle:
        for chunk in iter(lambda: file_handle.read(8192), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def group_replicate_runs(sample_ids: list[str]) -> dict[str, list[str]]:
    """Agrupa una lista de sample_id por muestra base, quedandose solo con
    los que siguen la convencion "{base}_run{n}". Las corridas de cada grupo
    quedan ordenadas por numero de corrida."""
    runs_by_base: dict[str, list[tuple[int, str]]] = {}
    for sample_id in sample_ids:
        parsed = parse_replicate_run_id(sample_id)
        if parsed is None:
            continue
        base_sample_id, run_number = parsed
        runs_by_base.setdefault(base_sample_id, []).append((run_number, sample_id))

    return {
        base_sample_id: [sample_id for _, sample_id in sorted(runs)]
        for base_sample_id, runs in runs_by_base.items()
    }


def build_gene_sets_by_run(amr_table: pd.DataFrame) -> dict[str, set[str]]:
    """Conjunto de genes detectados con confianza, por sample_id (corrida)."""
    confident_detections = amr_table.loc[amr_table["meets_identity_coverage_threshold"]]
    return confident_detections.groupby("sample_id")["gene_symbol"].apply(set).to_dict()


def build_assembly_hashes_by_run(run_ids: list[str], assemblies_dir: Path) -> dict[str, str | None]:
    """Hash SHA-256 del ensamblaje filtrado final de cada corrida. None si el
    archivo no existe (por ejemplo, si esa corrida todavia no se proceso)."""
    hashes = {}
    for run_id in run_ids:
        assembly_path = assemblies_dir / run_id / "contigs.filtered.fasta"
        hashes[run_id] = compute_sha256(assembly_path) if assembly_path.is_file() else None
    return hashes


def build_reproducibility_report(
    runs_by_base: dict[str, list[str]],
    gene_sets_by_run: dict[str, set[str]],
    status_by_run: dict[str, str],
    assemblies_dir: Path,
) -> pd.DataFrame:
    """Arma una fila por cada PAR de corridas dentro de cada muestra base,
    comparando genes, hash del ensamblaje y estado final."""
    comparison_rows = []

    for base_sample_id, run_ids in runs_by_base.items():
        if len(run_ids) < 2:
            continue  # no hay nada que comparar con una sola corrida

        assembly_hashes_by_run = build_assembly_hashes_by_run(run_ids, assemblies_dir)

        for run_a, run_b in itertools.combinations(run_ids, 2):
            genes_a = gene_sets_by_run.get(run_a, set())
            genes_b = gene_sets_by_run.get(run_b, set())
            hash_a = assembly_hashes_by_run.get(run_a)
            hash_b = assembly_hashes_by_run.get(run_b)
            status_a = status_by_run.get(run_a)
            status_b = status_by_run.get(run_b)

            comparison_rows.append({
                "base_sample_id": base_sample_id,
                "run_a": run_a,
                "run_b": run_b,
                "genes_run_a": ", ".join(sorted(genes_a)) if genes_a else "none",
                "genes_run_b": ", ".join(sorted(genes_b)) if genes_b else "none",
                "exact_gene_concordance": exact_gene_concordance(genes_a, genes_b),
                "jaccard_similarity": round(jaccard_similarity(genes_a, genes_b), 4),
                "identical_assembly_file": (
                    hash_a is not None and hash_a == hash_b
                ),
                "status_run_a": status_a,
                "status_run_b": status_b,
                "concordant_final_status": status_a is not None and status_a == status_b,
            })

    return pd.DataFrame(comparison_rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comparar corridas repetidas de la misma muestra para evaluar reproducibilidad."
    )
    parser.add_argument(
        "--amr-table", type=Path, default=Path("results/tables/amr_summary.tsv"),
        help="Listado largo de AMR de todas las corridas (results/tables/amr_summary.tsv)",
    )
    parser.add_argument(
        "--master-results", type=Path, default=Path("results/tables/master_results.tsv"),
        help="Tabla maestra de todas las corridas, con columna final_status",
    )
    parser.add_argument(
        "--assemblies-dir", type=Path, default=Path("results/assemblies"),
        help="Carpeta con los ensamblajes filtrados de cada corrida",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/tables/reproducibility_report.tsv"),
        help="Ruta de salida del reporte de reproducibilidad",
    )
    args = parser.parse_args()

    amr_table = pd.read_csv(args.amr_table, sep="\t")
    master_table = pd.read_csv(args.master_results, sep="\t")

    runs_by_base = group_replicate_runs(master_table["sample_id"].tolist())
    if not runs_by_base:
        raise ValueError(
            "No se encontraron sample_id con la convencion '{base}_runN' en "
            f"{args.master_results}; no hay corridas repetidas que comparar."
        )

    gene_sets_by_run = build_gene_sets_by_run(amr_table)
    status_by_run = master_table.set_index("sample_id")["final_status"].to_dict()

    reproducibility_report = build_reproducibility_report(
        runs_by_base, gene_sets_by_run, status_by_run, args.assemblies_dir
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    reproducibility_report.to_csv(args.output, sep="\t", index=False)

    fully_concordant = (
        (reproducibility_report["exact_gene_concordance"] == 1.0)
        & reproducibility_report["identical_assembly_file"]
        & reproducibility_report["concordant_final_status"]
    ).sum()
    print(
        f"{len(reproducibility_report)} par(es) de corridas comparados en "
        f"{len(runs_by_base)} muestra(s) base, escritos en {args.output}"
    )
    print(f"Pares totalmente concordantes (genes + archivo + estado): {fully_concordant}/{len(reproducibility_report)}")


if __name__ == "__main__":
    main()
