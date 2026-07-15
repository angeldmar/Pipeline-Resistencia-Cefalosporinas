# Pipeline de detección de resistencia a cefalosporinas de tercera generación en *E. coli*

## Descripción general

Este proyecto diseña, programa y evalúa un pipeline bioinformático de código abierto
para analizar secuencias genómicas de *Escherichia coli* asociadas con resistencia a
cefalosporinas de tercera generación. El flujo integra validación de datos, control
de calidad, ensamblaje genómico, identificación taxonómica, detección de genes de
resistencia y generación de reportes.

Su finalidad es demostrar **factibilidad técnica**, **desempeño computacional**,
**concordancia analítica** y **reproducibilidad** usando genomas obtenidos de
repositorios públicos. El proyecto **no pretende validación clínica** ni su
implementación inmediata en entornos institucionales: los reportes describen
determinantes genotípicos detectados, nunca conclusiones clínicas.

## Arquitectura

El pipeline sigue un diseño **Python-first**, con roles claramente separados:

- **Python**: validación de entradas y metadatos, descarga y organización de
  archivos, ejecución de herramientas bioinformáticas externas (vía `subprocess`),
  integración de resultados, generación de reportes, control de errores y
  trazabilidad (versiones, hashes, fechas).
- **Snakemake**: orquestación del flujo completo — cada regla declara sus
  entradas, salidas y comandos, y Snakemake construye las dependencias entre pasos.
- **R**: reservado únicamente para el análisis estadístico final (sensibilidad,
  especificidad, índice kappa, intervalos de confianza, coeficiente de variación,
  gráficas).

Todos los umbrales y parámetros que puedan cambiar (calidad, ensamblaje, taxonomía,
AMR, etc.) viven en archivos YAML dentro de `config/`, nunca fijados dentro de los
scripts.

## Estructura del proyecto

```
.
├── Snakefile
├── config/
│   ├── config.yaml               # parámetros centrales (umbrales, hilos, rutas)
│   ├── samples.tsv                # tabla de muestras (entrada principal del pipeline)
│   └── resistance_targets.yaml    # clasificación de familias de beta-lactamasas
├── workflow/
│   ├── rules/                     # reglas de Snakemake (.smk), una por etapa
│   ├── scripts/                   # scripts Python y R que hacen el trabajo real
│   └── envs/                      # un ambiente conda por herramienta externa
├── data/
│   ├── raw/                       # FASTQ crudos descargados
│   ├── reference/                 # bases de datos de referencia (Kraken2, AMRFinder, CheckM)
│   └── metadata/                  # manifiestos y metadatos generados
├── results/                       # salidas de cada etapa (qc, ensamblajes, taxonomía, amr, reportes, ...)
├── logs/                          # logs de ejecución de cada herramienta
└── tests/                         # unit/, integration/, e2e/, fixtures/
```

## Estado de avance

El desarrollo sigue un orden incremental: primero se programa y prueba manualmente
cada pieza en Python, y solo después se traslada a reglas de Snakemake con sus
ambientes Conda correspondientes (esto facilita ubicar el origen de cualquier error).

| # | Parte | Estado |
|---|-------|--------|
| 1 | Estructura de carpetas del proyecto | ✅ Hecho |
| 2 | Tabla de muestras (`config/samples.tsv`) | ✅ Hecho |
| 3 | Validación de metadatos (`validate_samples.py`) | ✅ Hecho |
| 4 | Configuración central (`config/config.yaml`) | ✅ Hecho |
| 5 | Clasificación de genes de resistencia (`config/resistance_targets.yaml`) | ✅ Hecho |
| 6 | Descarga de datos (`download_data.py`) | ✅ Hecho |
| 7 | Procesar manualmente una muestra real | ⏳ Pendiente |
| 8 | Control de calidad de lecturas (fastp) | ✅ Hecho |
| 9 | Verificación de cobertura | ✅ Hecho |
| 10 | Ensamblaje genómico (SPAdes) | ✅ Hecho |
| 11 | Evaluación del ensamblaje (QUAST) | ✅ Hecho |
| 12 | Completitud y contaminación (CheckM) | ✅ Hecho |
| 13 | Identificación taxonómica (Kraken2) | ✅ Hecho |
| 14 | Anotación genómica | ⏳ Pendiente |
| 15 | Detección de AMR (AMRFinderPlus) | ⏳ Pendiente |
| 16 | Comparación con estándar de referencia | ⏳ Pendiente |
| 17 | Integración de resultados (tabla maestra) | ⏳ Pendiente |
| 18 | Medición de desempeño computacional | ⏳ Pendiente |
| 19 | Pruebas de reproducibilidad | ⏳ Pendiente |
| 20 | Estadística en R | ⏳ Pendiente |
| 21 | Generación de reportes HTML | ⏳ Pendiente |
| 22 | Snakefile principal y reglas de Snakemake | ⏳ Pendiente |
| 23 | Ambientes Conda | ⏳ Pendiente |
| 24 | Pruebas (unitarias, integración, e2e, negativas) | ⏳ Pendiente |

## Detalle de lo implementado

### 1. Estructura de carpetas

Se creó la estructura fija del proyecto directamente en la raíz del repositorio
(ver árbol arriba). Las carpetas de datos/resultados que empiezan vacías incluyen
un `.gitkeep` para que git las versione. `tests/` se organizó en `unit/`,
`integration/`, `e2e/` y `fixtures/`, anticipando los 4 niveles de prueba
requeridos más adelante.

### 2. `config/samples.tsv`

Tabla de muestras con un identificador propio y estable (`sample_id`), independiente
del nombre de archivo. Columnas:

| Columna | Descripción |
|---|---|
| `sample_id` | identificador único de la muestra |
| `run_accession` | accesión SRA/ENA de la corrida de secuenciación |
| `biosample` | accesión BioSample |
| `sequencing_platform` | plataforma de secuenciación |
| `phenotype_cefotaxime` / `phenotype_ceftriaxone` / `phenotype_ceftazidime` | fenotipo (`S`/`I`/`R`/`NA`) por cefalosporina |
| `expected_genes` | gen o mecanismo de resistencia esperado (estándar de referencia) |
| `data_source` | fuente documentada del dato (trazabilidad) |

### 3. `workflow/scripts/validate_samples.py`

Valida `samples.tsv` **antes** de cualquier descarga o procesamiento, y junta
todos los errores encontrados en un solo reporte (en vez de detenerse en el
primero). Verifica:

- que no falte ninguna columna obligatoria;
- que no haya `sample_id` duplicados;
- que los fenotipos solo contengan `S`, `I`, `R` o `NA`;
- que `run_accession` y `biosample` tengan un formato de accesión reconocible
  (`SRR/ERR/DRR` + dígitos, `SAMN/SAMEA/SAMD` + dígitos);
- que cada muestra tenga una fuente de datos documentada (`data_source`);
- opcionalmente, si la tabla trae columnas `local_fastq_r1`/`local_fastq_r2`
  (muestras con lecturas ya locales en vez de descargarse), que esos archivos
  existan en disco.

Probado con la tabla válida (pasa) y con una tabla corrupta a propósito
(detecta duplicado, fenotipo inválido, accesiones mal formadas y fuente faltante).

### 4. `config/config.yaml`

Centraliza todos los parámetros configurables: hilos por herramienta, umbrales de
calidad/cobertura, umbrales de ensamblaje (QUAST), umbrales de completitud/
contaminación, umbrales taxonómicos (Kraken2), umbrales de AMR (AMRFinderPlus),
número de corridas de reproducibilidad y rutas a bases de datos de referencia.

### 5. `config/resistance_targets.yaml`

Clasifica familias de beta-lactamasas relevantes para cefalosporinas de tercera
generación en tres categorías: BLEE (`extended_spectrum_beta_lactamases`), AmpC
(`ampc`) y carbapenemasas (`carbapenemases`). Incluye una nota explícita de que el
prefijo de familia (p. ej. `TEM`, `SHV`) **no** basta para clasificar un gen como
BLEE — se necesita la subclase reportada por AMRFinderPlus o una tabla curada de
alelos.

### 6. `workflow/scripts/download_data.py`

Descarga las lecturas paired-end de cada muestra desde SRA usando `fasterq-dump`
(Python solo orquesta el proceso vía `subprocess`, no reemplaza la herramienta
oficial). Por cada muestra:

1. descarga con `fasterq-dump --split-files`;
2. verifica que se obtuvieron **ambos** archivos paired-end (falla si falta uno);
3. comprime a `data/raw/{sample_id}_R1.fastq.gz` / `_R2.fastq.gz`;
4. calcula el hash SHA-256 de cada archivo comprimido;
5. registra todo en `data/metadata/download_manifest.tsv` (una fila por archivo:
   muestra, accesión, repositorio, fuente, fecha, ruta, tamaño, hash y estado).

Si una muestra falla, se marca `FAILED` con el motivo del error y el script
continúa con las demás (una descarga fallida no debe detener el resto del lote).
Al reprocesar una muestra, sus filas antiguas en el manifiesto se reemplazan, no
se duplican. Soporta `--sample-id` para procesar una sola muestra.

Probado con `fasterq-dump` ausente (falla controladamente, queda registrado como
`FAILED`) y con un `fasterq-dump` simulado (ruta exitosa completa: descarga →
compresión → hash → manifiesto).

### 8. Control de calidad de lecturas (fastp)

**`workflow/rules/quality_control.smk`** define dos reglas de Snakemake:

- `fastp`: recorta adaptadores y filtra por calidad/longitud las lecturas crudas
  de `data/raw/`, usando los umbrales `quality.minimum_length` y
  `quality.minimum_phred` de `config.yaml` (no quedan fijos en la regla). Produce
  los FASTQ recortados en `results/trimmed/` y un reporte JSON + HTML en
  `results/qc/fastp/`.
- `parse_fastp`: llama a `parse_fastp.py parse` sobre el JSON de cada muestra.

**`workflow/scripts/parse_fastp.py`** reduce el JSON de fastp a las métricas que
le interesan al pipeline (lecturas iniciales/retenidas, % de lecturas filtradas,
bases Q20/Q30, % de contenido GC, % de duplicación). Tiene dos subcomandos:

- `parse <sample_id> <fastp.json>` → escribe una fila en
  `results/tables/fastp/{sample_id}.tsv` (un archivo **por muestra**).
- `combine` → junta todas las tablas por muestra en
  `results/tables/fastp_summary.tsv`.

**Decisión de arquitectura:** en vez de que cada muestra escriba directamente
sobre un único `fastp_summary.tsv` compartido (como sugiere literalmente la
sección 6.3 del diseño), cada muestra escribe su propio archivo. Snakemake puede
ejecutar varias muestras en paralelo, y varios procesos escribiendo a la vez
sobre el mismo archivo compartido causaría condiciones de carrera y resultados
corruptos o incompletos. La tabla combinada se arma en un paso de agregación
aparte (`combine`), que se conectará a la regla `all` del Snakefile cuando se
arme el flujo completo (parte de integración de resultados / Snakefile principal).
Este mismo patrón (archivo por muestra + agregación final) se reutilizará para
QUAST, taxonomía y AMRFinderPlus más adelante.

Probado con reportes JSON de fastp simulados (dos muestras con distintas tasas
de filtrado, duplicación y GC): las métricas extraídas y los porcentajes
calculados coinciden con el cálculo manual, y `combine` junta correctamente
ambas tablas individuales en un resumen ordenado por `sample_id`.

### 9. Verificación de cobertura

En vez de crear un script nuevo, la estimación de cobertura se integró
directamente en **`parse_fastp.py`**, porque el JSON de fastp ya trae todo lo
necesario para calcularla (lecturas retenidas y bases totales tras el
filtrado) — no hace falta ejecutar ninguna herramienta adicional.

- `estimate_coverage(read_count, mean_read_length, genome_size)`: implementa
  la fórmula `Cobertura = (lecturas × longitud media) / tamaño del genoma`,
  tal como la define el diseño del pipeline. La longitud media de lectura se
  calcula a partir de `total_bases / total_reads` (después del filtrado), en
  vez de asumirse fija.
- `classify_coverage(estimated_coverage, minimum_coverage, warning_coverage)`:
  clasifica la cobertura en `PASS` (≥30x), `WARNING` (15x–30x) o `FAIL` (<15x).
  Los umbrales 30x/15x están en `config.yaml`
  (`quality.minimum_coverage` / `quality.warning_coverage`); el documento solo
  especificaba el umbral PASS, así que el umbral WARNING (15x) se acordó
  explícitamente con el usuario antes de implementarlo.
- Ninguna muestra se descarta en silencio: `mean_read_length`,
  `estimated_coverage` y `coverage_status` quedan como columnas nuevas en la
  misma tabla por muestra (`results/tables/fastp/{sample_id}.tsv`), visibles
  también en el resumen combinado.

Probado con la fórmula del propio ejemplo del diseño (2,000,000 lecturas × 150
pb / 5,000,000 pb = 60x) y con tres reportes JSON simulados que producen
exactamente los tres estados (`PASS` a 60x, `WARNING` a 21x, `FAIL` a 9x).

### 10. Ensamblaje genómico (SPAdes)

**`workflow/rules/assembly.smk`** define dos reglas:

- `spades`: ensambla de novo las lecturas ya recortadas (`results/trimmed/`)
  usando el modo `--careful` (reduce errores de ensamblaje a costa de más
  tiempo de cómputo), produciendo `results/assemblies/{sample}/contigs.fasta`.
- `filter_contigs`: descarta los contigs más cortos que
  `assembly.minimum_contig_length` (config.yaml, 500 pb por defecto), dejando
  `results/assemblies/{sample}/contigs.filtered.fasta` — el archivo que usarán
  QUAST y AMRFinderPlus más adelante.

**`workflow/scripts/filter_contigs.py`** (script nuevo, no listado explícitamente
en la estructura fija de carpetas, pero requerido por la sección 8.3 del
diseño) implementa `filter_contigs()` con Biopython: lee el FASTA, conserva
solo los contigs con longitud ≥ al mínimo configurado, y avisa (sin fallar)
si un ensamblaje se queda sin contigs válidos, para que quede visible en los
logs en vez de fallar silenciosamente más adelante en QUAST.

Probado con un FASTA sintético de 5 contigs (mezcla de longitudes, incluyendo
un caso límite de exactamente 500 pb): retiene correctamente los 3 contigs
≥500 pb y descarta los 2 cortos. También probado el caso extremo de un
ensamblaje sin ningún contig válido: no falla, solo emite la advertencia.

### 11. Evaluación del ensamblaje (QUAST)

Se agregaron dos reglas más a **`workflow/rules/assembly.smk`** (QUAST evalúa
directamente la salida de `filter_contigs`, así que se mantiene en el mismo
archivo en vez de crear uno nuevo):

- `quast`: corre QUAST sobre `contigs.filtered.fasta` y genera
  `results/qc/quast/{sample}/report.tsv` (junto con el resto de los reportes
  de QUAST en esa misma carpeta).
- `parse_quast`: llama a `parse_quast.py parse` sobre ese `report.tsv`.

**`workflow/scripts/parse_quast.py`** sigue el mismo patrón que
`parse_fastp.py` (subcomandos `parse`/`combine`, archivo por muestra para
evitar condiciones de carrera en Snakemake). Extrae de QUAST: número de
contigs, contig más largo, longitud total, % GC y N50, y clasifica el
ensamblaje con `classify_assembly()`:

- **FAIL** si hay más contigs que `assembly.maximum_contigs`, o si la
  longitud total cae fuera del rango `[minimum_total_length,
  maximum_total_length]` esperado para *E. coli* — en ambos casos sin
  importar qué tan bueno sea el N50.
- **WARNING** si el N50 es menor a `assembly.n50_warning_threshold`
  (ensamblaje más fragmentado de lo ideal, pero utilizable).
- **PASS** en cualquier otro caso.

Todos estos umbrales vienen de `config.yaml` (ya centralizados desde la parte
4), no están fijos en el script.

**Nota técnica:** QUAST repite varias métricas con distintos umbrales de
longitud en el mismo `report.tsv` (ej. `# contigs (>= 1000 bp)` junto a
`# contigs` a secas). El parser usa coincidencia exacta de nombre de fila
(`# contigs`, `Total length`, `N50`, etc.) para no confundir la métrica
global del ensamblaje con las variantes filtradas por longitud.

Probado con tres `report.tsv` sintéticos (formato real de QUAST) que
representan los tres estados: uno con pocos contigs y N50 alto (PASS), uno
con N50 bajo pero todo lo demás normal (WARNING), y uno con exceso de
contigs a pesar de tener un N50 razonable, confirmando que el criterio FAIL
tiene prioridad sobre WARNING.

### 12. Completitud y contaminación (CheckM)

El documento dejaba la herramienta sin especificar ("una herramienta
específica"); se acordó con el usuario usar **CheckM**, el estándar de facto
para completitud/contaminación de genomas bacterianos vía genes marcadores de
copia única específicos del linaje. Se agregaron dos reglas más a
`workflow/rules/assembly.smk` (evalúa el mismo ensamblaje filtrado que QUAST):

- `checkm`: CheckM espera una **carpeta** de bins (genomas), no un archivo
  suelto, así que la regla copia primero el ensamblaje filtrado de la muestra
  a su propia carpeta de bins y luego corre `checkm lineage_wf --tab_table`.
- `parse_checkm`: llama a `parse_checkm.py parse` sobre la tabla resultante.

**`workflow/scripts/parse_checkm.py`** (script nuevo, no listado en la
estructura fija) sigue el mismo patrón de archivo-por-muestra +
`combine`, con una diferencia importante: el subcomando `combine` no solo
escribe `checkm_summary.tsv`, sino que además genera
**`results/tables/checkm_exclusions.tsv`**, un registro aparte con únicamente
las muestras `FAIL` y el motivo — tal como pide la sección 10 ("las muestras
que fallen pueden excluirse del análisis principal, pero deben permanecer en
el registro de exclusiones"). Ninguna muestra desaparece del resumen general.

Clasificación: `PASS` solo si completitud ≥ `assembly.minimum_completeness`
(95%) **y** contaminación < `assembly.maximum_contamination` (5%); cualquier
otro caso es `FAIL` (el documento no pide un nivel `WARNING` intermedio aquí).

Añadido también `workflow/envs/checkm.yaml` (placeholder, se completará en la
parte de ambientes Conda) y `threads.checkm` en `config.yaml`, ya que CheckM
tampoco estaba en la lista original de ambientes.

Probado con tres reportes CheckM sintéticos (formato `--tab_table` real):
un caso PASS, un caso que falla solo por completitud baja, y un caso que
falla solo por contaminación alta — ambos modos de fallo se detectan por
separado correctamente, y el registro de exclusiones lista ambos con su
motivo sin que desaparezcan del resumen combinado.

### 13. Identificación taxonómica (Kraken2)

**`workflow/rules/taxonomy.smk`** (archivo nuevo) define dos reglas:

- `kraken2`: clasifica taxonómicamente las lecturas ya recortadas
  (`results/trimmed/`) contra la base de datos de referencia
  (`paths.kraken_database`). Se corre sobre lecturas, no sobre el ensamblaje,
  que es el uso estándar de Kraken2 y evita que un ensamblaje ya colapsado
  enmascare una mezcla de especies.
- `parse_kraken2`: llama a `parse_kraken2.py parse` sobre el reporte.

**`workflow/scripts/parse_kraken2.py`** (script nuevo) extrae: taxón
predominante, % asignado a *E. coli*, y % de otras especies (contaminación),
y clasifica con `classify_taxonomy()` — misma lógica del documento (PASS si
≥90% *E. coli* y <5% contaminantes; WARNING si ≥70% *E. coli*; FAIL en el
resto), con los umbrales en `config.yaml` (`taxonomy.*`, ya centralizados
desde la parte 4).

**Regla de revisión manual para *Shigella*** (sección 11, requisito explícito
del documento): *Shigella* es genómicamente tan cercana a *E. coli* que
históricamente se consideraron la misma especie. Si Kraken2 asigna lecturas a
*Shigella*, el script **no** las suma al porcentaje de "otras especies" que
dispararía un FAIL — en vez de eso, marca `requires_manual_review=True` y dejo
un registro aparte (`results/tables/taxonomy_manual_review.tsv`) con esas
muestras, para que un analista decida caso por caso. Esto implementa
literalmente la instrucción del documento de no excluir automáticamente estos
aislamientos sin revisión.

Probado con tres reportes Kraken2 sintéticos (formato real, jerárquico):

- 95% *E. coli* sin contaminantes → **PASS**.
- 85% *E. coli* + 10% *Shigella* → **WARNING** (85% no alcanza el umbral PASS
  de 90%, independientemente de Shigella) y queda marcada para revisión
  manual — confirma que la regla de Shigella no "rescata" artificialmente el
  puntaje, solo evita que se cuente como contaminación.
- 40% *E. coli* con 45% *Klebsiella pneumoniae* → **FAIL**, y el taxón
  predominante se identifica correctamente como *Klebsiella pneumoniae* (no
  se asume que el taxón predominante sea siempre *E. coli*).

## Próximos pasos

Continuar con la **parte 14**: anotación genómica (sección 12 del diseño del
pipeline).
