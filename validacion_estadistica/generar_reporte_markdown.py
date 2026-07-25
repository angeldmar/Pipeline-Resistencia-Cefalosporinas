"""Genera la version Markdown del "Analisis de resultados" de la validacion.

Produce validacion_estadistica/reporte/analisis_resultados.md a partir de los
mismos artefactos que ya calculo el notebook (resultados/tables y
resultados/estadisticas), y copia las figuras necesarias a
reporte/figuras/ para que el documento se vea completo directamente en
GitHub (con imagenes y enlaces relativos). Es la contraparte navegable en el
repositorio del reporte en Word (generar_reporte_tesis.py); comparte su
contenido y no recalcula estadistica.

Uso:
    python generar_reporte_markdown.py
Salida:
    reporte/analisis_resultados.md
    reporte/figuras/*.png
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

VALIDATION_DIR = Path(__file__).resolve().parent
TABLES_DIR = VALIDATION_DIR / "resultados" / "tables"
STATS_DIR = VALIDATION_DIR / "resultados" / "estadisticas"
PLOTS_DIR = STATS_DIR / "plots"
SAMPLES_CSV = VALIDATION_DIR / "muestras" / "muestras_validacion_completo.csv"
REPORT_DIR = VALIDATION_DIR / "reporte"
FIGURES_DIR = REPORT_DIR / "figuras"
OUTPUT_PATH = REPORT_DIR / "analisis_resultados.md"

FIGURES = [
    "confusion_matrix_reporte.png",
    "metricas_comparadas.png",
    "cv_reproducidabilidad_reporte.png",  # nombre corregido abajo si no existe
]


def pct(x: float) -> str:
    return f"{x * 100:.1f}%"


def df_to_markdown(dataframe: pd.DataFrame, headers: list[str]) -> str:
    """Tabla Markdown (GitHub-flavored) desde un DataFrame, sin dependencias."""
    lines = ["| " + " | ".join(headers) + " |",
             "| " + " | ".join("---" for _ in headers) + " |"]
    for _, row in dataframe.iterrows():
        cells = ["" if pd.isna(value) else str(value) for value in row]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def copy_figures() -> dict[str, str]:
    """Copia a reporte/figuras las figuras del reporte y devuelve el mapeo
    nombre_logico -> ruta relativa (o None si la figura no existe)."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    wanted = {
        "confusion": "confusion_matrix_reporte.png",
        "metricas": "metricas_comparadas.png",
        "cv": "cv_reproducidabilidad_reporte.png",
    }
    # el generador del reporte guarda el CV como cv_reproducidabilidad_reporte.png
    # (ver generar_reporte_tesis.py); si no esta, se intenta el nombre del notebook.
    resolved = {}
    for key, filename in wanted.items():
        source = PLOTS_DIR / filename
        if key == "cv" and not source.is_file():
            source = PLOTS_DIR / "cv_reproducibilidad.png"
        if source.is_file():
            shutil.copy(source, FIGURES_DIR / source.name)
            resolved[key] = f"figuras/{source.name}"
        else:
            resolved[key] = None
    return resolved


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    figures = copy_figures()

    summary = pd.read_csv(TABLES_DIR / "validation_summary.tsv", sep="\t")
    meta = pd.read_csv(SAMPLES_CSV)
    mc = pd.read_csv(STATS_DIR / "classification_metrics.csv").set_index("metric")
    ma = pd.read_csv(STATS_DIR / "classification_metrics_all_samples.csv").set_index("metric")
    kappa = pd.read_csv(STATS_DIR / "kappa.csv").iloc[0]
    cv_path = STATS_DIR / "cv_reproducibilidad.csv"
    cv = pd.read_csv(cv_path) if cv_path.exists() else None

    n_total = len(summary)
    n_otra = int((summary["species_check"] == "otra_especie").sum())
    n_conf = int(mc.loc["Accuracy", "n"])
    n_all = int(ma.loc["Accuracy", "n"])

    out: list[str] = []
    out.append("# Análisis de resultados — Validación estadística")
    out.append("")
    out.append("> Documento generado a partir de los resultados del pipeline sobre el "
               "conjunto de validación. Es la versión navegable del reporte en Word "
               "(`resultados/reporte/`). Los números provienen de los artefactos que "
               "calcula `notebooks/validacion_estadistica.ipynb`; este documento no "
               "recalcula estadística.")
    out.append("")
    out.append(f"La validación contrastó la predicción genotípica del pipeline con el "
               f"estándar de referencia documentado para {n_total} genomas de *E. coli* "
               "obtenidos de repositorios públicos, repartidos entre controles positivos, "
               "negativos y límite (Tabla 1). El análisis examina primero la identidad "
               "real de las muestras, luego la concordancia entre predicción y "
               "referencia, y por último el origen de cada discrepancia.")
    out.append("")

    # Tabla 1
    composition = pd.crosstab(meta["Clasificación"], meta["Tipo de control"]).reset_index()
    composition.columns = ["Categoría"] + list(composition.columns[1:])
    out.append("**Tabla 1.** Composición del conjunto de validación por categoría de "
               "resistencia y tipo de control.")
    out.append("")
    out.append(df_to_markdown(composition, list(composition.columns)))
    out.append("")

    # Identidad de especie
    out.append("## Identidad de especie")
    out.append("")
    out.append(f"La verificación taxonómica con Kraken2 cambió la lectura de todo el "
               f"conjunto. **Catorce de los {n_total} genomas ({pct(n_otra/n_total)}) no "
               "eran *E. coli***, pese a figurar como tales en los metadatos de origen "
               "(Tabla 2). Los taxones detectados no se limitaban a parientes próximos "
               "como *Klebsiella* o *Salmonella*: aparecieron géneros tan alejados como "
               "*Listeria* y *Staphylococcus*, e incluso una secuencia de origen humano. "
               "Ninguna correspondía a los controles negativos ni a los límite; todas se "
               "habían incorporado como positivos. Su presencia hundía la sensibilidad "
               "aparente sin que mediara ningún fallo de detección, porque el gen buscado "
               "no existía en el organismo secuenciado. El análisis principal las "
               "descarta y conserva, en paralelo, el cálculo sin exclusión para dejar "
               "expuesto su efecto.")
    out.append("")
    otra = summary[summary["species_check"] == "otra_especie"][
        ["sample_id", "category", "confusion_category", "predominant_taxon"]
    ]
    out.append(f"**Tabla 2.** Las {n_otra} muestras cuya especie predominante no es "
               "*E. coli* según Kraken2.")
    out.append("")
    out.append(df_to_markdown(otra, ["Accesión", "Categoría", "Clasificación", "Taxón predominante (Kraken2)"]))
    out.append("")

    # Concordancia
    out.append("## Concordancia entre genotipo y referencia")
    out.append("")
    out.append(f"Sobre las {n_conf} muestras de especie confirmada con referencia "
               "evaluable, la matriz de confusión concentra los aciertos en la diagonal: "
               "41 verdaderos positivos y 26 verdaderos negativos, frente a 4 falsos "
               "positivos y 2 falsos negativos (Figura 1).")
    out.append("")
    if figures["confusion"]:
        out.append(f"![Matriz de confusión]({figures['confusion']})")
        out.append("")
        out.append(f"*__Figura 1.__ Matriz de confusión del pipeline frente al estándar "
                   f"de referencia sobre las {n_conf} muestras de especie confirmada "
                   "(VP: verdadero positivo; VN: verdadero negativo; FP: falso positivo; "
                   "FN: falso negativo).*")
        out.append("")
    out.append(f"La sensibilidad alcanzó el **{pct(mc.loc['Sensitivity','estimate'])}** "
               f"(IC 95%: {pct(mc.loc['Sensitivity','ci_lower'])}–{pct(mc.loc['Sensitivity','ci_upper'])}) "
               f"y la especificidad el **{pct(mc.loc['Specificity','estimate'])}** "
               f"(IC 95%: {pct(mc.loc['Specificity','ci_lower'])}–{pct(mc.loc['Specificity','ci_upper'])}), "
               f"con una exactitud del **{pct(mc.loc['Accuracy','estimate'])}** (Tabla 3). "
               "Los intervalos se calcularon por el método exacto de Clopper-Pearson. El "
               "contraste con el conjunto sin depurar es nítido: la sensibilidad caía al "
               f"{pct(ma.loc['Sensitivity','estimate'])} y la exactitud al "
               f"{pct(ma.loc['Accuracy','estimate'])} (Figura 2). Esa brecha no mide un "
               "cambio en el pipeline, sino la retirada del sesgo que introducían las "
               "muestras mal identificadas.")
    out.append("")
    metrics_table = pd.DataFrame({
        "Métrica": ["Sensibilidad", "Especificidad", "Exactitud"],
        f"Especie confirmada (n={n_conf})": [
            f"{pct(mc.loc[m,'estimate'])} ({pct(mc.loc[m,'ci_lower'])}–{pct(mc.loc[m,'ci_upper'])})"
            for m in ["Sensitivity", "Specificity", "Accuracy"]
        ],
        f"Todas las evaluables (n={n_all})": [
            f"{pct(ma.loc[m,'estimate'])} ({pct(ma.loc[m,'ci_lower'])}–{pct(ma.loc[m,'ci_upper'])})"
            for m in ["Sensitivity", "Specificity", "Accuracy"]
        ],
    })
    out.append("**Tabla 3.** Métricas de desempeño con intervalo de confianza del 95% "
               "(Clopper-Pearson), con y sin exclusión taxonómica.")
    out.append("")
    out.append(df_to_markdown(metrics_table, list(metrics_table.columns)))
    out.append("")
    if figures["metricas"]:
        out.append(f"![Métricas comparadas]({figures['metricas']})")
        out.append("")
        out.append("*__Figura 2.__ Métricas de desempeño antes y después de excluir las "
                   "muestras de otra especie. Las barras de error representan el IC del 95%.*")
        out.append("")
    out.append(f"El índice kappa de Cohen fija la concordancia corregida por azar en "
               f"**{kappa['kappa_value']:.3f} (p < 0.001)**, dentro del rango casi "
               "perfecto de la escala de Landis y Koch. Sin la depuración taxonómica "
               "descendía a 0.596, apenas moderado. La distancia entre ambos valores mide "
               "el peso de los catorce genomas espurios sobre la concordancia global.")
    out.append("")

    # Casos discordantes
    out.append("## Casos discordantes")
    out.append("")
    out.append("Los seis desacuerdos restantes admiten lectura individual (Tabla 4). Los "
               "cuatro falsos positivos aparecen todos en cepas *E. coli* negativas y "
               "comparten patrón: dos portan `blaTEM-1`, una β-lactamasa de espectro "
               "estrecho, y dos presentan mutaciones puntuales en `cirA`. Ninguno "
               "confiere por sí solo resistencia a cefalosporinas de tercera generación. "
               "El desacuerdo nace de la definición del control —ausencia de β-lactamasas "
               "de espectro extendido o AmpC adquiridas— frente a una regla de "
               "comparación que señala cualquier β-lactamasa. Es una discrepancia de "
               "criterio, no de detección. Los dos falsos negativos genuinos, confirmados "
               "como *E. coli* por MLST, no mostraron rastro del gen esperado ni por "
               "debajo de los umbrales de identidad y cobertura, y se reservan para "
               "revisión manual.")
    out.append("")
    disc = summary[summary["confusion_category"].isin(["FN", "FP"])][
        ["sample_id", "category", "expected_gene", "detected_gene", "confusion_category", "species_check"]
    ].copy()
    disc["species_check"] = disc["species_check"].map(
        {"confirmada": "E. coli", "otra_especie": "otra especie", "shigella_revision_manual": "Shigella"}
    )
    out.append("**Tabla 4.** Casos discordantes con su verificación de especie.")
    out.append("")
    out.append(df_to_markdown(disc, ["Accesión", "Categoría", "Gen esperado", "Gen detectado", "Clasif.", "Especie"]))
    out.append("")

    # Controles limite
    out.append("## Desempeño ante ensamblajes de calidad variable")
    out.append("")
    out.append("Los cinco controles límite, escogidos por la calidad de su ensamblaje y "
               "no por su genotipo, delimitan el rango operativo del método. La cepa del "
               "brote alemán de 2011 (O104:H4 TY-2482), con el ensamblaje más fragmentado "
               "del conjunto, todavía recuperó `blaCTX-M` y `blaTEM-1`, aunque perdió la "
               "resolución del alelo exacto y del tipo de secuencia. En el extremo "
               "opuesto, la cepa de referencia EC958 —clon pandémico ST131, con "
               "ensamblaje completo— entregó los cuatro genes descritos en la literatura. "
               "Las tres muestras intermedias, una de ellas la cepa tipo de la especie, "
               "no arrojaron genes adquiridos, y ambos motores de detección coincidieron "
               "en ese resultado. El desempeño decae de forma gradual con la calidad del "
               "ensamblaje, sin quiebres abruptos.")
    out.append("")

    # Reproducibilidad
    out.append("## Reproducibilidad computacional")
    out.append("")
    if cv is not None:
        cv_time_mean = cv["cv_tiempo_pct"].mean()
        cv_time_max = cv["cv_tiempo_pct"].max()
        cv_ram_mean = cv["cv_ram_pct"].mean()
        cv_ram_max = cv["cv_ram_pct"].max()
        out.append(f"Cada una de las cinco muestras de control se ejecutó tres veces "
                   f"(Tabla 5). El tiempo de ejecución se mantuvo estable entre corridas, "
                   f"con un coeficiente de variación promedio del {cv_time_mean:.1f}% y un "
                   f"máximo del {cv_time_max:.1f}%. La RAM máxima fluctuó más —CV promedio "
                   f"del {cv_ram_mean:.1f}%, hasta el {cv_ram_max:.1f}% en una muestra—, "
                   "porque el pico de memoria depende del paso de colocación filogenética "
                   "de CheckM, sensible a la presión de memoria del sistema y a la carga "
                   "concurrente del equipo durante cada corrida (Figura 3). El tiempo de "
                   "cómputo del pipeline es reproducible; el consumo de memoria pico "
                   "admite una variación moderada según las condiciones de ejecución.")
        out.append("")
        cv_table = cv.copy()
        cv_table.columns = ["Muestra", "Tiempo medio (s)", "CV tiempo (%)", "RAM media (GB)", "CV RAM (%)"]
        out.append("**Tabla 5.** Media y coeficiente de variación del tiempo de ejecución "
                   "y de la RAM máxima entre las tres corridas de cada muestra.")
        out.append("")
        out.append(df_to_markdown(cv_table, list(cv_table.columns)))
        out.append("")
        if figures["cv"]:
            out.append(f"![Coeficiente de variación]({figures['cv']})")
            out.append("")
            out.append("*__Figura 3.__ Coeficiente de variación del tiempo de ejecución y "
                       "de la RAM máxima entre las tres corridas de cada muestra de control.*")
            out.append("")
    else:
        out.append("_Pendiente: el cálculo del coeficiente de variación aún no está "
                   "disponible._")
        out.append("")

    OUTPUT_PATH.write_text("\n".join(out), encoding="utf-8")
    print(f"Reporte Markdown generado en {OUTPUT_PATH}")
    print(f"Figuras copiadas a {FIGURES_DIR}")


if __name__ == "__main__":
    main()
