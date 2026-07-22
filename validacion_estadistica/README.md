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

Pendiente: correr el pipeline sobre las 92 muestras y completar el
notebook de análisis.
