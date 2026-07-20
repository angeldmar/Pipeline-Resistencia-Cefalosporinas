"""Generacion del reporte HTML individual de una muestra.

Junta todo lo que ya calcularon los modulos anteriores (calidad, cobertura,
ensamblaje, taxonomia, completitud/contaminacion, genes de resistencia,
comparacion con la referencia, desempeno, versiones de herramientas) en un
unico reporte HTML autocontenido (la grafica va embebida en base64, no como
archivo aparte).

REGLA DE ALCANCE: el reporte describe determinantes GENOTIPICOS detectados,
nunca una conclusion clinica. En vez de "Aislado resistente a ceftriaxona",
el reporte dice algo como "Se detecto
el determinante blaCTX-M-15, asociado con beta-lactamasas de espectro
extendido (BLEE)...". Ver build_gene_interpretation_sentences().
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import base64
import io

import matplotlib
matplotlib.use("Agg")  # backend sin interfaz grafica, necesario para correr en servidores/CI
import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Environment, FileSystemLoader

# Descripcion en lenguaje genotipico (no clinico) de cada categoria
# mecanistica de beta-lactamasa, usada para construir las oraciones de
# interpretacion. "Other" cubre familias no clasificadas en BLEE/AmpC/
# carbapenemasa por classify_cephalosporin_genes.py.
BETA_LACTAMASE_CATEGORY_DESCRIPTIONS = {
    "ESBL": "beta-lactamasas de espectro extendido (BLEE), asociadas con resistencia a "
            "cefalosporinas de tercera generación",
    "AmpC": "beta-lactamasas tipo AmpC, asociadas con resistencia a cefalosporinas",
    "Carbapenemase": "carbapenemasas, con actividad hidrolítica frente a carbapenémicos y otros "
                      "beta-lactámicos",
    "Other": "un mecanismo de resistencia a beta-lactámicos no clasificado en las categorías "
             "BLEE/AmpC/carbapenemasa",
}

# Mensajes de advertencia por estado de modulo, en collect_warnings().
QC_GATE_LABELS = {
    "coverage_status": "cobertura de secuenciación",
    "assembly_status": "calidad del ensamblaje",
    "completeness_status": "completitud/contaminación del ensamblaje",
    "taxonomy_status": "identificación taxonómica",
}


def clean_context_value(value):
    """Reemplaza valores nulos/NaN (de columnas que no se pudieron unir para
    esta muestra) por un texto explicito, en vez de que Jinja2 los muestre
    como "nan" en el HTML."""
    if pd.isna(value):
        return "N/D"
    return value


def load_master_row(master_results_path: Path, sample_id: str) -> dict:
    """Carga la fila de la tabla maestra correspondiente a una muestra."""
    master_table = pd.read_csv(master_results_path, sep="\t")
    matching_rows = master_table.loc[master_table["sample_id"] == sample_id]
    if matching_rows.empty:
        raise ValueError(f"No se encontro la muestra '{sample_id}' en {master_results_path}")
    raw_row = matching_rows.iloc[0].to_dict()
    return {key: clean_context_value(value) for key, value in raw_row.items()}


def load_amr_genes(amr_classified_path: Path, sample_id: str) -> list[dict]:
    """Carga las filas de genes de resistencia detectados para una muestra."""
    amr_table = pd.read_csv(amr_classified_path, sep="\t")
    sample_genes = amr_table.loc[amr_table["sample_id"] == sample_id]
    return sample_genes.to_dict("records")


def load_tool_versions(tool_versions_path: Path) -> dict:
    """Carga las versiones de herramientas (ver capture_tool_versions.py).
    Devuelve un diccionario vacio si el archivo todavia no existe."""
    if not tool_versions_path.is_file():
        return {}
    versions_table = pd.read_csv(tool_versions_path, sep="\t")
    return dict(zip(versions_table["tool"], versions_table["version"]))


def load_performance_by_module(performance_summary_path: Path, sample_id: str) -> list[dict]:
    """Carga el desglose de tiempo/CPU/RAM por modulo de una muestra."""
    if not performance_summary_path.is_file():
        return []
    performance_table = pd.read_csv(performance_summary_path, sep="\t")
    sample_rows = performance_table.loc[performance_table["sample_id"] == sample_id]
    return sample_rows.to_dict("records")


def build_gene_interpretation_sentences(amr_genes: list[dict]) -> list[str]:
    """Construye una oracion de interpretacion GENOTIPICA por cada gen
    detectado con confianza (nunca una conclusion clinica: ver docstring
    del modulo)."""
    sentences = []
    for gene in amr_genes:
        if not gene.get("meets_identity_coverage_threshold"):
            continue  # las detecciones de baja confianza no se interpretan

        gene_symbol = gene["gene_symbol"]
        percent_identity = gene["percent_identity"]
        percent_coverage = gene["percent_coverage"]

        if str(gene.get("antimicrobial_class", "")).upper() == "BETA-LACTAM":
            category = gene.get("beta_lactamase_category", "Other")
            category_description = BETA_LACTAMASE_CATEGORY_DESCRIPTIONS.get(
                category, BETA_LACTAMASE_CATEGORY_DESCRIPTIONS["Other"]
            )
            sentences.append(
                f"Se detectó el determinante {gene_symbol}, asociado con {category_description} "
                f"(identidad {percent_identity}%, cobertura {percent_coverage}%)."
            )
        else:
            antimicrobial_class = str(gene.get("antimicrobial_class", "desconocida")).lower()
            sentences.append(
                f"Se detectó el determinante {gene_symbol}, asociado con resistencia a "
                f"{antimicrobial_class} (identidad {percent_identity}%, cobertura {percent_coverage}%)."
            )
    return sentences


def collect_warnings(master_row: dict, amr_genes: list[dict]) -> list[str]:
    """Junta, en lenguaje llano, todas las advertencias relevantes para esta
    muestra: estados WARNING/FAIL de cualquier modulo de QC, presencia de
    Shigella (revision manual), detecciones de AMR de baja confianza, y un
    estado final distinto de PASS."""
    warnings = []

    for status_column, label in QC_GATE_LABELS.items():
        status_value = master_row.get(status_column)
        if status_value in ("WARNING", "FAIL"):
            warnings.append(f"Estado {status_value} en {label}.")

    if master_row.get("requires_manual_review") in (True, "True"):
        warnings.append(
            "Se detectaron lecturas asignadas a Shigella; se recomienda revisión manual por la "
            "cercanía genómica entre Shigella y E. coli."
        )

    low_confidence_genes = [
        gene["gene_symbol"] for gene in amr_genes if not gene.get("meets_identity_coverage_threshold")
    ]
    if low_confidence_genes:
        warnings.append(
            "Detecciones por debajo del umbral de identidad/cobertura configurado (no "
            f"interpretadas arriba): {', '.join(low_confidence_genes)}."
        )

    final_status = master_row.get("final_status")
    if final_status not in ("PASS", None):
        warnings.append(f"Estado final de la muestra: {final_status}.")

    return warnings


def render_amr_confidence_chart(amr_genes: list[dict]) -> str | None:
    """Genera un grafico de barras de % identidad / % cobertura por gen
    detectado, devuelto como PNG codificado en base64 (para embeber
    directamente en el HTML, sin archivos de imagen sueltos). None si la
    muestra no tiene genes detectados."""
    if not amr_genes:
        return None

    gene_labels = [gene["gene_symbol"] for gene in amr_genes]
    identity_values = [gene["percent_identity"] for gene in amr_genes]
    coverage_values = [gene["percent_coverage"] for gene in amr_genes]

    bar_positions = range(len(gene_labels))
    bar_width = 0.35

    figure, axes = plt.subplots(figsize=(max(4, len(gene_labels) * 1.3), 4))
    axes.bar([position - bar_width / 2 for position in bar_positions], identity_values, bar_width, label="% Identidad")
    axes.bar([position + bar_width / 2 for position in bar_positions], coverage_values, bar_width, label="% Cobertura")
    axes.set_xticks(list(bar_positions))
    axes.set_xticklabels(gene_labels, rotation=30, ha="right")
    axes.set_ylim(0, 105)
    axes.set_ylabel("%")
    axes.set_title("Identidad y cobertura de los genes detectados")
    axes.legend()
    figure.tight_layout()

    image_buffer = io.BytesIO()
    figure.savefig(image_buffer, format="png", dpi=100)
    plt.close(figure)
    image_buffer.seek(0)
    return base64.b64encode(image_buffer.read()).decode("ascii")


def build_csv_download_data_uri(rows: list[dict]) -> str | None:
    """Convierte una lista de filas (diccionarios) en un enlace de descarga
    CSV embebido como data URI en base64, para que el reporte siga siendo un
    unico archivo HTML autocontenido (mismo criterio ya usado para el
    grafico embebido). None si no hay filas que exportar."""
    if not rows:
        return None
    csv_buffer = io.StringIO()
    pd.DataFrame(rows).to_csv(csv_buffer, index=False)
    encoded_csv = base64.b64encode(csv_buffer.getvalue().encode("utf-8")).decode("ascii")
    return f"data:text/csv;charset=utf-8;base64,{encoded_csv}"


def render_report(context: dict, template_dir: Path, template_name: str) -> str:
    """Renderiza la plantilla Jinja2 del reporte con el contexto dado."""
    environment = Environment(loader=FileSystemLoader(str(template_dir)))
    template = environment.get_template(template_name)
    return template.render(**context)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generar el reporte HTML individual de una muestra."
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--master-results", type=Path, default=Path("results/tables/master_results.tsv"))
    parser.add_argument("--amr-classified", type=Path, default=Path("results/tables/amr_classified.tsv"))
    parser.add_argument("--performance-summary", type=Path, default=Path("results/tables/performance_summary.tsv"))
    parser.add_argument("--tool-versions", type=Path, default=Path("data/metadata/tool_versions.tsv"))
    parser.add_argument("--template-dir", type=Path, default=Path("workflow/templates"))
    parser.add_argument("--template-name", type=str, default="sample_report.html.j2")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    master_row = load_master_row(args.master_results, args.sample_id)
    amr_genes = load_amr_genes(args.amr_classified, args.sample_id)
    tool_versions = load_tool_versions(args.tool_versions)
    performance_by_module = load_performance_by_module(args.performance_summary, args.sample_id)

    context = {
        **master_row,
        "amr_genes": amr_genes,
        "gene_interpretations": build_gene_interpretation_sentences(amr_genes),
        "warnings": collect_warnings(master_row, amr_genes),
        "amr_chart_base64": render_amr_confidence_chart(amr_genes),
        "tool_versions": tool_versions,
        "performance_by_module": performance_by_module,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "amr_genes_csv_uri": build_csv_download_data_uri(amr_genes),
        "sample_summary_csv_uri": build_csv_download_data_uri([master_row]),
    }

    html_content = render_report(context, args.template_dir, args.template_name)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html_content, encoding="utf-8")
    print(f"Reporte de {args.sample_id} escrito en {args.output}")


if __name__ == "__main__":
    main()
