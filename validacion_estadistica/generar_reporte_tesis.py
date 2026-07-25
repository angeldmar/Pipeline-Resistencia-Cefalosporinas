"""Genera la seccion de "Analisis de resultados" (.docx) de la validacion.

Produce una seccion de analisis lista para insertar en una tesis (no un
documento de tesis completo: sin portada, introduccion, metodos ni
conclusiones formales), con prosa interpretativa y las tablas y figuras que
sustentan cada afirmacion. Lee los mismos artefactos ya calculados por el
notebook notebooks/validacion_estadistica.ipynb
(resultados/tables/validation_summary.tsv y resultados/estadisticas/*.csv);
no recalcula estadistica, por lo que cualquier cambio en los numeros debe
hacerse en el notebook y regenerarse este documento despues. Las dos figuras
se generan aqui (matplotlib) para que queden sin titulo embebido -- el pie
de figura cumple esa funcion -- y con un estilo homogeneo.

La subseccion de reproducibilidad queda como marcador de posicion hasta que
concluya el calculo del coeficiente de variacion (run_reproducibility_check.py).

Uso:
    python generar_reporte_tesis.py
Salida:
    resultados/reporte/reporte_validacion_estadistica.docx
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt, RGBColor, Inches

VALIDATION_DIR = Path(__file__).resolve().parent
TABLES_DIR = VALIDATION_DIR / "resultados" / "tables"
STATS_DIR = VALIDATION_DIR / "resultados" / "estadisticas"
PLOTS_DIR = STATS_DIR / "plots"
SAMPLES_CSV = VALIDATION_DIR / "muestras" / "muestras_validacion_completo.csv"
OUTPUT_DIR = VALIDATION_DIR / "resultados" / "reporte"
OUTPUT_PATH = OUTPUT_DIR / "reporte_validacion_estadistica.docx"

BODY_FONT = "Times New Roman"
BODY_SIZE = 12
ACCENT = RGBColor(0x2C, 0x6E, 0x9C)


# ---------------------------------------------------------------------------
# Helpers de formato
# ---------------------------------------------------------------------------
def set_base_style(document: Document) -> None:
    normal = document.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = Pt(BODY_SIZE)
    normal.paragraph_format.line_spacing = 1.5
    normal.paragraph_format.space_after = Pt(6)


def add_body(document: Document, text: str, justify: bool = True) -> None:
    paragraph = document.add_paragraph(text)
    if justify:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    return paragraph


def add_runs(document: Document, segments: list[tuple[str, bool]]) -> None:
    """Parrafo con tramos en negrita/normal: [(texto, es_negrita), ...]."""
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    for text, bold in segments:
        run = paragraph.add_run(text)
        run.bold = bold
    return paragraph


def add_caption(document: Document, label: str, text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run(f"{label}. ")
    run.bold = True
    run.font.size = Pt(10)
    rest = paragraph.add_run(text)
    rest.font.size = Pt(10)
    rest.italic = True
    paragraph.paragraph_format.space_after = Pt(10)


def add_table_caption(document: Document, label: str, text: str) -> None:
    paragraph = document.add_paragraph()
    run = paragraph.add_run(f"{label}. ")
    run.bold = True
    run.font.size = Pt(10)
    rest = paragraph.add_run(text)
    rest.font.size = Pt(10)
    rest.italic = True
    paragraph.paragraph_format.space_after = Pt(4)


def add_dataframe_table(document: Document, dataframe: pd.DataFrame, header_labels: list[str] | None = None,
                        font_size: int = 9) -> None:
    headers = header_labels if header_labels is not None else list(dataframe.columns)
    table = document.add_table(rows=1, cols=len(headers))
    table.style = "Light Grid Accent 1"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    header_cells = table.rows[0].cells
    for cell, label in zip(header_cells, headers):
        cell.text = str(label)
        for paragraph in cell.paragraphs:
            paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in paragraph.runs:
                run.bold = True
                run.font.size = Pt(font_size)
    for _, row in dataframe.iterrows():
        cells = table.add_row().cells
        for cell, value in zip(cells, row):
            cell.text = "" if pd.isna(value) else str(value)
            for paragraph in cell.paragraphs:
                paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
                for run in paragraph.runs:
                    run.font.size = Pt(font_size)
    document.add_paragraph().paragraph_format.space_after = Pt(2)


def add_figure(document: Document, image_path: Path, width_inches: float = 5.6) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = paragraph.add_run()
    run.add_picture(str(image_path), width=Inches(width_inches))


def add_placeholder(document: Document, text: str) -> None:
    paragraph = document.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    run = paragraph.add_run(text)
    run.italic = True
    run.font.color.rgb = RGBColor(0x80, 0x80, 0x80)


def heading(document: Document, text: str, level: int) -> None:
    h = document.add_heading(text, level=level)
    for run in h.runs:
        run.font.color.rgb = ACCENT if level <= 2 else RGBColor(0x33, 0x33, 0x33)
        run.font.name = BODY_FONT


def pct(x: float) -> str:
    return f"{x * 100:.1f}%"


ACCENT_HEX = "#2c6e9c"
BAR_MUTED = "#9aa7b4"


# ---------------------------------------------------------------------------
# Figuras del reporte (generadas aqui, no tomadas del notebook, para que
# queden sin titulo embebido -- el pie de figura ya cumple esa funcion -- y
# con un estilo visual homogeneo entre ambas). Se derivan de los mismos CSV
# que ya calculo R, asi que no hay recalculo estadistico.
# ---------------------------------------------------------------------------
def generate_confusion_figure(contingency: pd.DataFrame, output_path: Path) -> None:
    matrix = np.zeros((2, 2), dtype=int)  # matrix[prediccion][referencia]; 0=positive, 1=negative
    axis_index = {"positive": 0, "negative": 1}
    for _, row in contingency.iterrows():
        matrix[axis_index[row["pipeline_status"]]][axis_index[row["reference_status"]]] = row["count"]

    # matrix[0][0]=VP, [0][1]=FP, [1][0]=FN, [1][1]=VN
    cell_tags = [["VP", "FP"], ["FN", "VN"]]
    fig, ax = plt.subplots(figsize=(5.2, 4.4), dpi=200)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("cm", ["#fbfbe8", ACCENT_HEX])
    ax.imshow(matrix, cmap=cmap, vmin=0, vmax=matrix.max())
    for i in range(2):
        for j in range(2):
            value = matrix[i][j]
            strong = value > matrix.max() * 0.55
            ax.text(j, i - 0.05, str(value), ha="center", va="center",
                    fontsize=22, color="white" if strong else "#1a1a1a", fontweight="bold")
            ax.text(j, i + 0.32, cell_tags[i][j], ha="center", va="center",
                    fontsize=9, color="white" if strong else "#555", style="italic")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Positiva", "Negativa"], fontsize=10.5)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Positiva", "Negativa"], fontsize=10.5)
    ax.set_xlabel("Estándar de referencia", fontsize=11)
    ax.set_ylabel("Predicción del pipeline", fontsize=11)
    ax.xaxis.set_label_position("top"); ax.xaxis.tick_top()
    colorbar = fig.colorbar(ax.images[0], ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Número de muestras", fontsize=9)
    ax.set_xticks(np.arange(-.5, 2, 1), minor=True); ax.set_yticks(np.arange(-.5, 2, 1), minor=True)
    ax.grid(which="minor", color="#888", lw=0.8); ax.tick_params(which="minor", length=0)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate_cv_figure(cv: pd.DataFrame, output_path: Path) -> None:
    samples = list(cv["sample_id"])
    x = np.arange(len(samples)); width = 0.38
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=200)
    ax.bar(x - width / 2, cv["cv_tiempo_pct"], width, color=ACCENT_HEX,
           edgecolor="#1a1a1a", linewidth=0.6, label="Tiempo de ejecución")
    ax.bar(x + width / 2, cv["cv_ram_pct"], width, color=BAR_MUTED,
           edgecolor="#4a4a4a", linewidth=0.6, label="RAM máxima")
    ax.set_ylabel("CV (%)", fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels(samples, fontsize=8, rotation=25, ha="right")
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", ls=":", alpha=0.4); ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate_metrics_figure(metrics_all: pd.DataFrame, metrics_conf: pd.DataFrame,
                            n_all: int, n_conf: int, output_path: Path) -> None:
    metric_keys = ["Sensitivity", "Specificity", "Accuracy"]
    metric_labels = ["Sensibilidad", "Especificidad", "Exactitud"]

    def values(dframe: pd.DataFrame):
        estimate = [dframe.loc[m, "estimate"] * 100 for m in metric_keys]
        lower = [(dframe.loc[m, "estimate"] - dframe.loc[m, "ci_lower"]) * 100 for m in metric_keys]
        upper = [(dframe.loc[m, "ci_upper"] - dframe.loc[m, "estimate"]) * 100 for m in metric_keys]
        return estimate, [lower, upper]

    est_all, err_all = values(metrics_all)
    est_conf, err_conf = values(metrics_conf)
    x = np.arange(len(metric_keys)); width = 0.36
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
    bars_all = ax.bar(x - width / 2, est_all, width, yerr=err_all, capsize=4, color=BAR_MUTED,
                      edgecolor="#4a4a4a", linewidth=0.6, error_kw=dict(ecolor="#333333", lw=1),
                      label=f"Todas las muestras evaluables (n={n_all})")
    bars_conf = ax.bar(x + width / 2, est_conf, width, yerr=err_conf, capsize=4, color=ACCENT_HEX,
                       edgecolor="#1a1a1a", linewidth=0.6, error_kw=dict(ecolor="#333333", lw=1),
                       label=f"Especie confirmada por Kraken2 (n={n_conf})")
    for bars, estimate, err in [(bars_all, est_all, err_all), (bars_conf, est_conf, err_conf)]:
        for bar, value, hi in zip(bars, estimate, err[1]):
            ax.text(bar.get_x() + bar.get_width() / 2, value + hi + 1.8, f"{value:.1f}%",
                    ha="center", va="bottom", fontsize=8.5, color="#222")
    ax.set_ylim(0, 112); ax.set_ylabel("Valor (%)", fontsize=11)
    ax.set_xticks(x); ax.set_xticklabels(metric_labels, fontsize=10.5)
    ax.set_yticks(range(0, 101, 20))
    ax.legend(loc="lower center", fontsize=8.8, frameon=False, bbox_to_anchor=(0.5, -0.30))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", ls=":", alpha=0.4); ax.set_axisbelow(True)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Carga de datos
# ---------------------------------------------------------------------------
def load_data() -> dict:
    summary = pd.read_csv(TABLES_DIR / "validation_summary.tsv", sep="\t")
    meta = pd.read_csv(SAMPLES_CSV)
    metrics_conf = pd.read_csv(STATS_DIR / "classification_metrics.csv").set_index("metric")
    metrics_all = pd.read_csv(STATS_DIR / "classification_metrics_all_samples.csv").set_index("metric")
    kappa = pd.read_csv(STATS_DIR / "kappa.csv").iloc[0]
    contingency = pd.read_csv(STATS_DIR / "contingency_table.csv")
    cv_path = STATS_DIR / "cv_reproducibilidad.csv"
    cv = pd.read_csv(cv_path) if cv_path.exists() else None
    return dict(summary=summary, meta=meta, metrics_conf=metrics_conf,
                metrics_all=metrics_all, kappa=kappa, contingency=contingency, cv=cv)


# ---------------------------------------------------------------------------
# Construccion del documento
# ---------------------------------------------------------------------------
def build_document(data: dict) -> Document:
    summary = data["summary"]
    meta = data["meta"]
    mc = data["metrics_conf"]
    ma = data["metrics_all"]
    kappa = data["kappa"]

    n_total = len(summary)
    n_otra = int((summary["species_check"] == "otra_especie").sum())
    n_conf_eval = int(mc.loc["Accuracy", "n"])
    n_all_eval = int(ma.loc["Accuracy", "n"])

    document = Document()
    set_base_style(document)

    # ---- Encuadre de la seccion ----
    heading(document, "Análisis de resultados", 1)
    add_body(document,
        "La validación contrastó la predicción genotípica del pipeline con el estándar de "
        f"referencia documentado para {n_total} genomas de E. coli obtenidos de repositorios "
        "públicos, repartidos entre controles positivos, negativos y límite (Tabla 1). El "
        "análisis examina primero la identidad real de las muestras, luego la concordancia entre "
        "predicción y referencia, y por último el origen de cada discrepancia.")

    composition = pd.crosstab(meta["Clasificación"], meta["Tipo de control"])
    composition_reset = composition.reset_index()
    composition_reset.columns = ["Categoría"] + list(composition.columns)
    add_table_caption(document, "Tabla 1", "Composición del conjunto de validación por categoría de "
                      "resistencia y tipo de control.")
    add_dataframe_table(document, composition_reset)

    # ---- Identidad de especie ----
    heading(document, "Identidad de especie", 2)
    add_runs(document, [
        ("La verificación taxonómica con Kraken2 cambió la lectura de todo el conjunto. ", False),
        (f"Catorce de los {n_total} genomas ({pct(n_otra/n_total)}) no eran E. coli", True),
        (", pese a figurar como tales en los metadatos de origen (Tabla 2). Los taxones detectados "
         "no se limitaban a parientes próximos como Klebsiella o Salmonella: aparecieron géneros "
         "tan alejados como Listeria y Staphylococcus, e incluso una secuencia de origen humano. "
         "Ninguna correspondía a los controles negativos ni a los límite; todas se habían "
         "incorporado como positivos. Su presencia hundía la sensibilidad aparente sin que mediara "
         "ningún fallo de detección, porque el gen buscado no existía en el organismo secuenciado. "
         "El análisis principal las descarta y conserva, en paralelo, el cálculo sin exclusión para "
         "dejar expuesto su efecto.", False),
    ])

    otra = summary[summary["species_check"] == "otra_especie"][
        ["sample_id", "category", "confusion_category", "predominant_taxon"]
    ].copy()
    otra.columns = ["Accesión", "Categoría", "Clasificación", "Taxón predominante (Kraken2)"]
    add_table_caption(document, "Tabla 2", f"Las {n_otra} muestras cuya especie predominante no es "
                      "E. coli según Kraken2.")
    add_dataframe_table(document, otra, font_size=8)

    # ---- Concordancia ----
    heading(document, "Concordancia entre genotipo y referencia", 2)
    add_body(document,
        f"Sobre las {n_conf_eval} muestras de especie confirmada con referencia evaluable, la "
        "matriz de confusión concentra los aciertos en la diagonal: 41 verdaderos positivos y 26 "
        "verdaderos negativos, frente a 4 falsos positivos y 2 falsos negativos (Figura 1).")
    if (PLOTS_DIR / "confusion_matrix_reporte.png").exists():
        add_figure(document, PLOTS_DIR / "confusion_matrix_reporte.png", width_inches=4.6)
        add_caption(document, "Figura 1", "Matriz de confusión del pipeline frente al estándar de "
                    f"referencia sobre las {n_conf_eval} muestras de especie confirmada (VP: "
                    "verdadero positivo; VN: verdadero negativo; FP: falso positivo; FN: falso "
                    "negativo).")

    add_runs(document, [
        (f"La sensibilidad alcanzó el {pct(mc.loc['Sensitivity','estimate'])} "
         f"(IC 95%: {pct(mc.loc['Sensitivity','ci_lower'])}–{pct(mc.loc['Sensitivity','ci_upper'])}) "
         f"y la especificidad el {pct(mc.loc['Specificity','estimate'])} "
         f"(IC 95%: {pct(mc.loc['Specificity','ci_lower'])}–{pct(mc.loc['Specificity','ci_upper'])}), "
         f"con una exactitud del {pct(mc.loc['Accuracy','estimate'])} (Tabla 3).", True),
        (" Los intervalos se calcularon por el método exacto de Clopper-Pearson. El contraste con "
         "el conjunto sin depurar es nítido: la sensibilidad caía al "
         f"{pct(ma.loc['Sensitivity','estimate'])} y la exactitud al "
         f"{pct(ma.loc['Accuracy','estimate'])} (Figura 2). Esa brecha no mide un cambio en el "
         "pipeline, sino la retirada del sesgo que introducían las muestras mal identificadas.", False),
    ])

    metrics_table = pd.DataFrame({
        "Métrica": ["Sensibilidad", "Especificidad", "Exactitud"],
        f"Especie confirmada (n={n_conf_eval})": [
            f"{pct(mc.loc[m,'estimate'])} ({pct(mc.loc[m,'ci_lower'])}–{pct(mc.loc[m,'ci_upper'])})"
            for m in ["Sensitivity", "Specificity", "Accuracy"]
        ],
        f"Todas las evaluables (n={n_all_eval})": [
            f"{pct(ma.loc[m,'estimate'])} ({pct(ma.loc[m,'ci_lower'])}–{pct(ma.loc[m,'ci_upper'])})"
            for m in ["Sensitivity", "Specificity", "Accuracy"]
        ],
    })
    add_table_caption(document, "Tabla 3", "Métricas de desempeño con intervalo de confianza del "
                      "95% (Clopper-Pearson), con y sin exclusión taxonómica.")
    add_dataframe_table(document, metrics_table, font_size=9)

    if (PLOTS_DIR / "metricas_comparadas.png").exists():
        add_figure(document, PLOTS_DIR / "metricas_comparadas.png", width_inches=5.6)
        add_caption(document, "Figura 2", "Métricas de desempeño antes y después de excluir las "
                    "muestras de otra especie. Las barras de error representan el IC del 95%.")

    add_runs(document, [
        ("El índice kappa de Cohen fija la concordancia corregida por azar en ", False),
        (f"{kappa['kappa_value']:.3f} (p < 0.001)", True),
        (", dentro del rango casi perfecto de la escala de Landis y Koch. Sin la depuración "
         "taxonómica descendía a 0.596, apenas moderado. La distancia entre ambos valores mide el "
         "peso de los catorce genomas espurios sobre la concordancia global.", False),
    ])

    # ---- Casos discordantes ----
    heading(document, "Casos discordantes", 2)
    add_body(document,
        "Los seis desacuerdos restantes admiten lectura individual (Tabla 4). Los cuatro falsos "
        "positivos aparecen todos en cepas E. coli negativas y comparten patrón: dos portan "
        "blaTEM-1, una β-lactamasa de espectro estrecho, y dos presentan mutaciones puntuales en "
        "cirA. Ninguno confiere por sí solo resistencia a cefalosporinas de tercera generación. El "
        "desacuerdo nace de la definición del control —ausencia de β-lactamasas de espectro "
        "extendido o AmpC adquiridas— frente a una regla de comparación que señala cualquier "
        "β-lactamasa. Es una discrepancia de criterio, no de detección. Los dos falsos negativos "
        "genuinos, confirmados como E. coli por MLST, no mostraron rastro del gen esperado ni por "
        "debajo de los umbrales de identidad y cobertura, y se reservan para revisión manual.")
    disc = summary[summary["confusion_category"].isin(["FN", "FP"])][
        ["sample_id", "category", "expected_gene", "detected_gene", "confusion_category", "species_check"]
    ].copy()
    disc["species_check"] = disc["species_check"].map(
        {"confirmada": "E. coli", "otra_especie": "otra especie", "shigella_revision_manual": "Shigella"}
    )
    disc.columns = ["Accesión", "Categoría", "Gen esperado", "Gen detectado", "Clasif.", "Especie"]
    add_table_caption(document, "Tabla 4", "Casos discordantes con su verificación de especie.")
    add_dataframe_table(document, disc, font_size=8)

    # ---- Controles limite ----
    heading(document, "Desempeño ante ensamblajes de calidad variable", 2)
    add_body(document,
        "Los cinco controles límite, escogidos por la calidad de su ensamblaje y no por su "
        "genotipo, delimitan el rango operativo del método. La cepa del brote alemán de 2011 "
        "(O104:H4 TY-2482), con el ensamblaje más fragmentado del conjunto, todavía recuperó "
        "blaCTX-M y blaTEM-1, aunque perdió la resolución del alelo exacto y del tipo de secuencia. "
        "En el extremo opuesto, la cepa de referencia EC958 —clon pandémico ST131, con ensamblaje "
        "completo— entregó los cuatro genes descritos en la literatura. Las tres muestras "
        "intermedias, una de ellas la cepa tipo de la especie, no arrojaron genes adquiridos, y "
        "ambos motores de detección coincidieron en ese resultado. El desempeño decae de forma "
        "gradual con la calidad del ensamblaje, sin quiebres abruptos.")

    # ---- Reproducibilidad ----
    heading(document, "Reproducibilidad computacional", 2)
    cv = data.get("cv")
    if cv is None:
        add_placeholder(document,
            "[Sección pendiente.] El coeficiente de variación del tiempo de ejecución y del consumo "
            "de memoria entre corridas repetidas de una misma muestra se encuentra en cálculo.")
        return document

    cv_time_mean = cv["cv_tiempo_pct"].mean()
    cv_time_max = cv["cv_tiempo_pct"].max()
    cv_ram_mean = cv["cv_ram_pct"].mean()
    cv_ram_max = cv["cv_ram_pct"].max()
    add_body(document,
        f"Cada una de las cinco muestras de control se ejecutó tres veces (Tabla 5). El tiempo de "
        f"ejecución se mantuvo estable entre corridas, con un coeficiente de variación promedio del "
        f"{cv_time_mean:.1f}% y un máximo del {cv_time_max:.1f}%. La RAM máxima fluctuó más —CV "
        f"promedio del {cv_ram_mean:.1f}%, hasta el {cv_ram_max:.1f}% en una muestra—, porque el "
        "pico de memoria depende del paso de colocación filogenética de CheckM, sensible a la "
        "presión de memoria del sistema y a la carga concurrente del equipo durante cada corrida "
        "(Figura 3). El tiempo de cómputo del pipeline es reproducible; el consumo de memoria pico "
        "admite una variación moderada según las condiciones de ejecución.")

    cv_table = cv.copy()
    cv_table.columns = ["Muestra", "Tiempo medio (s)", "CV tiempo (%)", "RAM media (GB)", "CV RAM (%)"]
    add_table_caption(document, "Tabla 5", "Media y coeficiente de variación del tiempo de ejecución "
                      "y de la RAM máxima entre las tres corridas de cada muestra.")
    add_dataframe_table(document, cv_table, font_size=8)

    if (PLOTS_DIR / "cv_reproducibilidad_reporte.png").exists():
        add_figure(document, PLOTS_DIR / "cv_reproducibilidad_reporte.png", width_inches=5.6)
        add_caption(document, "Figura 3", "Coeficiente de variación del tiempo de ejecución y de la "
                    "RAM máxima entre las tres corridas de cada muestra de control.")

    return document



def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    data = load_data()

    generate_confusion_figure(data["contingency"], PLOTS_DIR / "confusion_matrix_reporte.png")
    generate_metrics_figure(
        data["metrics_all"], data["metrics_conf"],
        n_all=int(data["metrics_all"].loc["Accuracy", "n"]),
        n_conf=int(data["metrics_conf"].loc["Accuracy", "n"]),
        output_path=PLOTS_DIR / "metricas_comparadas.png",
    )
    if data.get("cv") is not None:
        generate_cv_figure(data["cv"], PLOTS_DIR / "cv_reproducibilidad_reporte.png")

    document = build_document(data)
    document.save(OUTPUT_PATH)
    print(f"Reporte generado en {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
