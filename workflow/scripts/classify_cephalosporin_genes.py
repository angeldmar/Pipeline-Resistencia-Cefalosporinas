"""Clasificacion mecanistica de beta-lactamasas detectadas por AMRFinderPlus.

Toma el listado largo ya normalizado (parse_amrfinder.py) y, para cada gen
de la clase BETA-LACTAM, decide a que categoria mecanistica pertenece segun
config/resistance_targets.yaml: BLEE (beta-lactamasa de espectro extendido),
AmpC, carbapenemasa, u "Other" si no coincide con ninguna familia conocida.

A diferencia del ejemplo original del diseno del pipeline (que fijaba las
familias directamente en el codigo), esta clasificacion LEE las familias
desde resistance_targets.yaml, para poder agregar o quitar familias sin
tocar este script. Tambien generaliza la regla de revision de subclase: no
solo SHV/TEM necesitan que la subclase confirme "EXTENDED-SPECTRUM" antes de
llamarlos BLEE (ver nota en resistance_targets.yaml sobre por que el solo
prefijo del simbolo no basta) -- GES tiene el mismo problema (existen alelos
GES que son carbapenemasas, no BLEE), asi que cualquier familia marcada con
el sufijo "-ESBL" en el YAML pasa por la misma verificacion de subclase.
"""

from __future__ import annotations

from pathlib import Path
import argparse

import pandas as pd
import yaml

# Sufijo usado en resistance_targets.yaml para marcar una familia cuyo
# prefijo NO basta por si solo para llamarla BLEE: se necesita ademas que la
# subclase reportada por AMRFinderPlus confirme "EXTENDED-SPECTRUM".
SUBCLASS_CONFIRMATION_SUFFIX = "-ESBL"
EXTENDED_SPECTRUM_SUBCLASS_MARKER = "EXTENDED-SPECTRUM"

BETA_LACTAM_CLASS_NAME = "BETA-LACTAM"


def load_resistance_targets(resistance_targets_path: Path) -> dict:
    """Carga las familias de beta-lactamasas clasificadas por categoria."""
    with open(resistance_targets_path) as targets_file:
        return yaml.safe_load(targets_file)


def classify_beta_lactamase(gene_symbol: str, antimicrobial_subclass: str, resistance_targets: dict) -> str:
    """Clasifica un gen de beta-lactamasa en 'ESBL', 'AmpC', 'Carbapenemase'
    u 'Other', usando las familias definidas en resistance_targets.yaml.

    El simple prefijo de familia (ej. "SHV") no basta para las familias
    marcadas con "-ESBL" en el YAML: se exige que la subclase reportada por
    AMRFinderPlus confirme "EXTENDED-SPECTRUM" antes de llamarlas BLEE, ya
    que la misma familia tiene alelos de espectro estrecho (o, en el caso de
    GES, incluso carbapenemasas).
    """
    symbol_upper = gene_symbol.upper()
    subclass_upper = (antimicrobial_subclass or "").upper()

    for esbl_family in resistance_targets.get("extended_spectrum_beta_lactamases", []):
        requires_subclass_confirmation = esbl_family.endswith(SUBCLASS_CONFIRMATION_SUFFIX)
        family_prefix = (
            esbl_family[: -len(SUBCLASS_CONFIRMATION_SUFFIX)]
            if requires_subclass_confirmation
            else esbl_family
        )
        if family_prefix.upper() not in symbol_upper:
            continue
        if not requires_subclass_confirmation:
            return "ESBL"
        if EXTENDED_SPECTRUM_SUBCLASS_MARKER in subclass_upper:
            return "ESBL"
        # El prefijo coincide (ej. "SHV") pero la subclase no confirma que
        # este alelo especifico sea de espectro extendido: se sigue
        # evaluando el resto de categorias en vez de asumir "Other" de una.

    for carbapenemase_family in resistance_targets.get("carbapenemases", []):
        if carbapenemase_family.upper() in symbol_upper:
            return "Carbapenemase"

    for ampc_family in resistance_targets.get("ampc", []):
        if ampc_family.upper() in symbol_upper:
            return "AmpC"

    return "Other"


def classify_amr_table(amr_table: pd.DataFrame, resistance_targets: dict) -> pd.DataFrame:
    """Agrega la columna beta_lactamase_category al listado largo de AMR.

    Solo se clasifican las filas de clase BETA-LACTAM; el resto (ej. genes de
    resistencia a aminoglucosidos) quedan marcadas como "N/A" porque la
    clasificacion BLEE/AmpC/carbapenemasa no les aplica.
    """
    classified_table = amr_table.copy()

    def classify_row(row: pd.Series) -> str:
        if str(row["antimicrobial_class"]).upper() != BETA_LACTAM_CLASS_NAME:
            return "N/A"
        return classify_beta_lactamase(row["gene_symbol"], row["antimicrobial_subclass"], resistance_targets)

    classified_table["beta_lactamase_category"] = classified_table.apply(classify_row, axis=1)
    return classified_table


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clasificar los genes de beta-lactamasa de un listado de AMR en BLEE/AmpC/carbapenemasa."
    )
    parser.add_argument(
        "amr_summary", type=Path,
        help="Listado largo de AMR ya normalizado (ej. results/tables/amr_summary.tsv)",
    )
    parser.add_argument(
        "--resistance-targets", type=Path, default=Path("config/resistance_targets.yaml"),
        help="Ruta a resistance_targets.yaml",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/tables/amr_classified.tsv"),
        help="Ruta de salida del listado clasificado",
    )
    args = parser.parse_args()

    amr_table = pd.read_csv(args.amr_summary, sep="\t")
    resistance_targets = load_resistance_targets(args.resistance_targets)
    classified_table = classify_amr_table(amr_table, resistance_targets)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    classified_table.to_csv(args.output, sep="\t", index=False)
    print(f"{len(classified_table)} deteccion(es) clasificadas, escritas en {args.output}")


if __name__ == "__main__":
    main()
