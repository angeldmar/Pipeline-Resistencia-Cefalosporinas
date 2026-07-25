# Validación estadística

Segunda sección del proyecto, distinta de la herramienta operativa del
pipeline (`Snakefile`, `workflow/`, `webapp/`). Mientras esa primera
sección se enfoca en correr el pipeline sobre una muestra dada y producir
su reporte, esta sección analiza, en conjunto y con muestras reales
aportadas para este fin, qué tan bien concuerdan las predicciones
genotípicas del pipeline con los fenotipos de referencia documentados.

## Estructura

- `muestras/`: metadatos, accesiones y genomas ensamblados (`genomas/`) de
  las muestras reales usadas para la validación. A diferencia de
  `data/raw/` en la sección operativa, estas SÍ se versionan: es un
  conjunto curado y acotado para esta validación, no un lote reproducible
  bajo demanda. Solo se ignora el Excel de accesiones de entrada (dato de
  entrada, no parte del conjunto ya validado).
- `resultados/`: resultados del pipeline recolectados para cada una de
  esas muestras (tablas, reportes), como entrada al análisis. Estos sí se
  regeneran corriendo el pipeline, mismo criterio que `results/` en la
  sección operativa.
- `notebooks/validacion_estadistica.ipynb`: notebook de R (kernel IRkernel)
  con tres apartados — muestras reales, resultados por muestra, y una
  discusión final que integra todos los casos.
- `generar_reporte_tesis.py`: genera un reporte en Word (formato académico,
  con sus tablas y figuras) a partir de los mismos artefactos que ya calculó
  el notebook. Es solo la capa de presentación para tesis; no recalcula
  estadística. La salida va a `resultados/reporte/` (no versionada, se
  regenera corriendo el script).

## Estado

**Conjunto de muestras definido y descargado: 92 genomas** de *E. coli*
reales, ensamblados, en `muestras/genomas/`. Criterio original: 100 genomas
(30 CTX-M, 15 SHV/TEM, 15 AmpC plasmídica, 30 negativos, 5 límite, 5 de
reproducibilidad tomados de los anteriores). Composición final lograda:

| Categoría | Objetivo | Logrado |
|---|---|---|
| CTX-M (ESBL) | 30 | 35 (30 + 5 de diversidad adicional para compensar el hueco de AmpC) |
| SHV/TEM (ESBL) | 15 | 12 |
| AmpC (CMY/DHA/FOX/ACC-MOX) | 15 | 10 |
| Negativos | 30 | 30 |
| Controles límite | 5 | 5 |
| **Total distinto** | 95 (máximo aritmético de los criterios) | **92** |

**Brechas y por qué no se completaron:** para 3 posiciones SHV/TEM y 1
AmpC, el biosample/accesión que traía el Excel original resultó inválido
(no existe en NCBI) o correspondía a un organismo/experimento
completamente distinto (un caso: RNA-seq humano etiquetado como *E. coli*
blaTEM-12). Se buscaron sustitutos reales en NCBI para varias de estas
posiciones, pero el pool de genomas de *E. coli* con blaSHV/blaTEM/AmpC
plasmídica **ya ensamblados** en NCBI es pequeño y se agotó rápido (las
búsquedas nuevas devolvían los mismos candidatos ya usados). Ante esto se
descartaron las posiciones irresolubles en vez de forzar sustitutos de
menor calidad — ver `Observaciones` en `muestras_validacion_completo.csv`
para el detalle fila por fila de cada sustitución y su justificación.

**Limitación de verificación de los negativos:** los 30 controles
negativos se verificaron por identidad de cepa (cepas de referencia
extensamente publicadas y caracterizadas: K-12, CFT073, O157:H7 Sakai,
Nissle 1917, BW25113, ATCC 25922, etc.), no volviendo a correr
AMRFinderPlus desde cero sobre cada una. La verificación genotípica real
con la versión de AMRFinderPlus de este proyecto ocurre naturalmente al
correr el pipeline sobre ellas — que es, de hecho, el propósito de esta
validación.

**Pipeline corrido, verificación taxonómica corrida, notebook completo y
ejecutado.** `run_validation_batch.py` procesó las 92 muestras (modo
solo-ensamblaje: QUAST, CheckM, AMRFinderPlus, ABricate, MLST, comparación
contra el estándar de referencia), y `run_taxonomy_check.py` corrió
Kraken2 sobre cada ensamblaje para confirmar la especie de cada muestra.

Ese segundo paso surgió al revisar los primeros falsos negativos: varios
tenían %GC muy alejado del rango normal de *E. coli* y MLST sin ningún
alelo resuelto. Kraken2 confirmó la sospecha — **14 de las 92 muestras
(15.2%) no son realmente *E. coli*** (Klebsiella, Salmonella, Citrobacter,
Listeria, Staphylococcus, Aeromonas, y una predominantemente humana),
casi con seguridad por sustituciones de accesión sin verificar durante el
armado del conjunto (ver `Brechas y por qué no se completaron` arriba).
Ninguna de las 14 es un control negativo ni un control límite.

`prepare_statistics_input.py` junta reference_comparison, CheckM, QUAST,
MLST, Kraken2 (`species_check`) y la categoría de control en
`resultados/tables/validation_summary.tsv`, y
`notebooks/validacion_estadistica.ipynb` (R vía IRkernel, ambiente
`workflow/envs/r_statistics.yaml`) calcula sensibilidad, especificidad,
exactitud (con IC 95%) y kappa de Cohen dos veces — sin filtrar por especie
(87 muestras evaluables) y excluyendo las 14 de otra especie (73 muestras)
— y desglosa los casos discordantes uno por uno. Sobre las 73 de especie
confirmada: **sensibilidad 95.3%** (IC 95% 84.2–99.4%), **especificidad
86.7%** (IC 95% 69.3–96.2%), **exactitud 91.8%** (IC 95% 83.0–96.9%),
**kappa 0.829** (concordancia casi perfecta) — frente a sensibilidad 77.2%
y kappa 0.596 sin ese filtro. Ver la sección "Discusión" del notebook para
el detalle completo.

`run_reproducibility_check.py` cierra la parte de reproducibilidad
computacional: corre 3 veces cada una de las 5 muestras marcadas como
control de reproducibilidad (midiendo tiempo y RAM por corrida, cada réplica
en su propio proceso para una medición limpia) y el notebook calcula el
coeficiente de variación. El tiempo de ejecución resultó estable (CV
promedio 4.0%) y la RAM máxima algo más variable (CV promedio 10.8%), por la
sensibilidad del pico de memoria de CheckM a la carga del sistema.

`generar_reporte_tesis.py` produce, a partir de los mismos artefactos, la
sección de "Análisis de resultados" en Word (`resultados/reporte/`), con su
prosa, tablas y figuras. El notebook también se puede exportar a PDF con
`jupyter nbconvert --to pdf` (requiere una distribución LaTeX, p. ej.
TinyTeX).
