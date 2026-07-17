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
| 14 | Anotación genómica (Prokka) | ✅ Hecho |
| 15 | Detección de AMR (AMRFinderPlus) | ✅ Hecho |
| 16 | Comparación con estándar de referencia | ✅ Hecho |
| 17 | Integración de resultados (tabla maestra) | ✅ Hecho |
| 18 | Medición de desempeño computacional | ✅ Hecho |
| 19 | Pruebas de reproducibilidad | ✅ Hecho |
| 20 | Estadística en R | ✅ Hecho |
| 21 | Generación de reportes HTML | ✅ Hecho |
| 22 | Snakefile principal y reglas de Snakemake | ✅ Hecho |
| 23 | Ambientes Conda | ✅ Hecho |
| 24 | Pruebas (unitarias, integración, extremo a extremo, negativas) | ✅ Hecho |

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
  (`quality.minimum_coverage` / `quality.warning_coverage`); la especificación
  del pipeline solo definía el umbral PASS, así que el umbral WARNING (15x) se
  fijó explícitamente antes de implementarlo.
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

La especificación del pipeline dejaba la herramienta sin definir ("una
herramienta específica"); se utiliza **CheckM**, el estándar de facto para
completitud/contaminación de genomas bacterianos vía genes marcadores de
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

### 14. Anotación genómica (Prokka)

Herramienta no especificada en la especificación del pipeline; se utiliza
**Prokka** (estándar clásico, ampliamente citado, sin necesidad de descargar
una base de datos de referencia grande, a diferencia de Bakta). Se creó
**`workflow/rules/annotation.smk`** (archivo nuevo, no estaba en la lista fija
de reglas — igual que ocurrió con algunos scripts en partes anteriores) con
dos reglas:

- `prokka`: anota el ensamblaje filtrado y produce `.gff`, `.gbk`, `.faa`,
  `.ffn`, además de `.tsv` y `.txt` (que Prokka genera igual, y que el
  parser usa para no tener que recalcular nada). Se le indica género/especie
  (`--genus Escherichia --species coli`) para una anotación más precisa.
- `parse_prokka`: llama a `parse_prokka.py parse` sobre el resumen de Prokka.

**`workflow/scripts/parse_prokka.py`** (script nuevo) **no vuelve a anotar
nada**: solo organiza lo que Prokka ya calculó. Extrae CDS, ARN ribosómico y
ARN de transferencia directamente del resumen `{sample}.txt` que Prokka
genera (evita recalcular conteos que la propia herramienta ya entrega),
cuenta los genes anotados como `"hypothetical protein"` (sin función
conocida) a partir de la tabla `{sample}.tsv`, y registra la versión de
Prokka consultando `prokka --version` en tiempo de ejecución (si la
herramienta no está disponible, registra `"unknown"` en vez de fallar, ya
que es un dato de trazabilidad, no bloqueante).

Se agregó `results/annotation/` a la estructura de resultados (no estaba en
la lista original de subcarpetas de `results/`) y `threads.prokka` en
`config.yaml`.

Como nota el propio documento, este paso es informativo/complementario:
AMRFinderPlus (parte 15) trabaja directamente sobre el ensamblaje de
nucleótidos y no depende de esta anotación.

Probado con un resumen y una tabla de Prokka sintéticos (formato real):
extrae correctamente CDS/rRNA/tRNA del resumen, cuenta bien 2 genes
hipotéticos de 5 filas de ejemplo, y maneja sin fallar el caso de Prokka no
instalado (versión `"unknown"`).

### 15. Detección de AMR (AMRFinderPlus)

**`workflow/rules/amr_detection.smk`** define tres reglas:

- `amrfinder`: corre AMRFinderPlus directamente sobre el ensamblaje de
  nucleótidos filtrado (no depende de la anotación de Prokka).
- `parse_amrfinder`: normaliza la salida cruda por muestra.
- `classify_cephalosporin_genes`: **paso de agregación final** (no por
  muestra, a diferencia de los demás parsers) — toma el listado ya combinado
  de todas las muestras y le agrega la clasificación mecanicista.

**`workflow/scripts/parse_amrfinder.py`** convierte la tabla cruda de
AMRFinderPlus a un esquema propio y estable (tabla **larga**: una fila por
gen detectado, siguiendo la recomendación de la sección 15 de evitar
combinar múltiples genes en columnas difíciles de analizar): símbolo del
gen, familia (derivada quitando el sufijo de alelo, ej. `blaCTX-M-15` →
`blaCTX-M`), alelo, nombre descriptivo, clase/subclase de antimicrobiano,
método de detección, % identidad, % cobertura, contig y coordenadas. Las
detecciones por debajo de los umbrales `amr.minimum_identity` /
`amr.minimum_gene_coverage` **no se descartan**: quedan marcadas con
`meets_identity_coverage_threshold=False`, visible en la tabla. Si una
muestra no tiene genes detectados, se escribe igual un archivo (con
encabezado, 0 filas) para distinguir "sin genes" de "no procesada".

**`workflow/scripts/classify_cephalosporin_genes.py`**: implementa
`classify_beta_lactamase()`, pero **generaliza** el ejemplo del documento en
un punto importante — en vez de fijar las familias BLEE/AmpC en el código
(duplicando lo que ya vive en `resistance_targets.yaml`), las **lee** de ese
archivo. Esto además corrigió dos huecos del ejemplo original: (1) GES tiene
el mismo problema que SHV/TEM (hay alelos GES que son carbapenemasas, no
BLEE), así que cualquier familia marcada `-ESBL` en el YAML pasa por la
misma verificación de subclase, no solo SHV y TEM explícitamente; (2) se
agregó la categoría `Carbapenemase` (el ejemplo original no la contemplaba
pese a que `resistance_targets.yaml` ya definía esa lista).

Probado con AMRFinderPlus sintético cubriendo los casos clave:

- `blaCTX-M-15` (subclase `CEPHALOSPORIN`) → **ESBL** — familia sin ambigüedad.
- `blaTEM-1` y `blaSHV-1` (subclase `BETA-LACTAM`, no `EXTENDED-SPECTRUM`) →
  **Other**, no ESBL — confirma la regla central de la sección 13.3: el
  prefijo de familia por sí solo no basta.
- `blaCMY-2` (subclase `CEPHALOSPORIN`) → **AmpC**.
- `blaNDM-1` (subclase `CARBAPENEM`) → **Carbapenemase**.
- Un gen de aminoglucósido → `N/A` (la clasificación BLEE no le aplica).
- Una detección con 60% de cobertura (< umbral 80%) → queda en la tabla con
  `meets_identity_coverage_threshold=False`, no se pierde.
- Una muestra sin genes detectados → tabla vacía sin error.

### 16. Comparación con el estándar de referencia

**`workflow/scripts/compare_to_reference.py`** (script nuevo) arma, por
muestra, la matriz TP/TN/FP/FN que R usará después para calcular
sensibilidad/especificidad/kappa — **ningún cálculo estadístico se hace
aquí**, solo se prepara y clasifica la matriz de datos, tal como pide el
documento. Compara `expected_genes` de `samples.tsv` contra los genes
detectados con confianza (`meets_identity_coverage_threshold=True`) en
`amr_summary.tsv`:

- Coincidencia **exacta** de alelo → `TP`.
- Coincidencia de **familia** pero alelo distinto al de referencia (ej. se
  esperaba `blaCTX-M-15` y se detectó `blaCTX-M-27`) → sigue contando como
  detección positiva (`TP`), pero queda registrado en `match_type=family`
  para que se note que el alelo exacto no coincidió.
- `expected_genes="none"` sin genes beta-lactámicos detectados → `TN`; si
  aparece alguno inesperado → `FP`.
- Gen esperado no detectado → `FN`.
- **`Indeterminado`** (quinta categoría que pide la sección 14, más allá del
  ejemplo de `confusion_category()` del documento): cuando `expected_genes`
  está vacío o `NA` — es decir, no hay estándar de referencia documentado
  para esa muestra — en vez de forzar una comparación sin sentido.

**Alcance deliberadamente limitado por ahora:** esta comparación es
únicamente a nivel de gen. Todavía no considera si la muestra falló algún
control de calidad anterior (cobertura, ensamblaje, taxonomía,
completitud) — esa integración le corresponde a la tabla maestra
(`merge_results.py`, parte 17), que es la que tiene visibilidad de todos los
módulos a la vez. Esta limitación de alcance queda documentada explícitamente
aquí, para la integración posterior.

**Bug encontrado y corregido durante las pruebas:** la primera versión trataba
`reference_positive` con `if reference_positive:`, y como `None` (el caso
indeterminado) también evalúa como falso en Python, cualquier muestra sin
referencia caía silenciosamente en la rama de "no se esperaba ningún gen" en
vez de su propia rama. Se corrigió comparando explícitamente los tres
estados (`is None` / verdadero / falso). De paso también apareció el mismo
problema de `pandas` leyendo el texto `"NA"` como valor nulo real (ya
resuelto antes en `validate_samples.py` con `.fillna("NA")`, pero no estaba
aplicado en este script nuevo) — corregido de la misma forma.

**Nota sobre el grafo de Snakemake:** las reglas `classify_cephalosporin_genes`
y `compare_to_reference` (agregadas a `amr_detection.smk`) dependen de
`results/tables/amr_summary.tsv`, que se arma juntando la tabla de *todas*
las muestras (`parse_amrfinder.py combine`). Esa regla de "combinar todas las
muestras" necesita la lista completa de muestras (`SAMPLES`), que recién
queda definida al construir el Snakefile principal — así que, por ahora,
`amr_summary.tsv` (y los demás `*_summary.tsv` de los módulos anteriores) se
generan manualmente vía CLI, y se conectarán al grafo completo cuando se
arme el Snakefile principal.

Probado: coincidencia exacta (TP), coincidencia por familia con alelo
distinto (TP, `match_type=family`), gen esperado no detectado (FN), ausencia
de gen confirmada (TN), y referencia indeterminada con detección real
(`Indeterminado`, con el gen detectado igual visible en la tabla en vez de
perderse).

### 17. Integración de resultados (tabla maestra)

**`workflow/scripts/merge_results.py`** arma `results/tables/master_results.tsv`:
**una fila por muestra**, uniendo (left join sobre `sample_id`, partiendo
siempre de `samples.tsv` como base para que ninguna muestra documentada
desaparezca) los resúmenes de todos los módulos: metadatos/fenotipo, calidad
y cobertura (fastp), métricas de ensamblaje (QUAST), completitud/contaminación
(CheckM), taxonomía (Kraken2), un resumen corto de genes de AMR detectados, y
la comparación con el estándar de referencia.

Siguiendo la advertencia de la sección 15 de evitar "combinar múltiples genes
en columnas difíciles de analizar", la tabla maestra **no** lista el detalle
de cada gen en columnas — solo agrega un resumen corto por muestra
(`detected_gene_count`, `detected_beta_lactam_gene_count`, `detected_genes`
como texto separado por comas). El detalle completo (identidad, cobertura,
coordenadas) sigue viviendo exclusivamente en la tabla larga
(`amr_summary.tsv` / `amr_classified.tsv`).

**`final_status`** (nueva columna, resuelve la limitación de alcance señalada
en la parte 16): combina los cuatro estados de control de calidad
(`coverage_status`, `assembly_status`, `completeness_status`,
`taxonomy_status`) en un solo veredicto — `EXCLUDED` si algún módulo dio
`FAIL` (candidata a excluirse del análisis principal, pero **sigue apareciendo
en la tabla**, nunca desaparece en silencio), `WARNING` si nada falló pero
algo advirtió, `PASS` si todo pasó, o `PENDING` si a esa muestra todavía no
le corrió ningún módulo de QC.

Columnas con el mismo nombre en dos módulos (`gc_content_percent` lo reporta
tanto fastp —de las lecturas— como QUAST —del ensamblaje—) se renombran antes
de unir (`reads_gc_content_percent` / `assembly_gc_content_percent`) para que
no se pisen entre sí.

**Columnas pendientes:** "Tiempo" y "Memoria" (medición de desempeño
computacional) aún no existen como módulo — se agregan en la parte 18,
siguiendo el orden recomendado del propio documento (primero la tabla
maestra, después tiempo/RAM). `merge_results.py` se extenderá entonces.

**Resiliencia:** si la tabla de resumen de algún módulo todavía no existe
(corrida parcial durante el desarrollo), se avisa por stdout y se continúa
sin ese módulo, en vez de fallar — probado quitando la tabla de CheckM: el
script avisa, sigue, y `final_status` se recalcula usando solo los módulos
efectivamente disponibles (una muestra que antes era `EXCLUDED` por
CheckM vuelve a `PASS` si ya no hay evidencia de ese fallo disponible).

Probado de extremo a extremo con tres muestras sintéticas consistentes entre
todos los módulos: una que pasa todo (`PASS`), una que falla completitud vía
CheckM (`EXCLUDED`), y una con cobertura límite (`WARNING`) — el resumen de
AMR y la comparación de referencia coinciden correctamente en cada fila.

### 18. Medición de desempeño computacional

**Decisión de arquitectura:** el documento sugiere envolver cada comando con
`/usr/bin/time -v` (GNU time). Esa variante `-v` es específica de GNU y no
existe en macOS/BSD (el entorno de desarrollo de este proyecto), ni se puede
asumir instalada en cualquier entorno de ejecución futuro. En su lugar se usó
el módulo `resource` de la biblioteca estándar de Python, que mide lo mismo
(tiempo, CPU, RAM máxima) de forma portable entre Linux y macOS, siguiendo el
mismo principio ya usado en todo el pipeline: Python orquesta las
herramientas externas, no depende de utilidades de shell específicas de una
plataforma.

**`workflow/scripts/run_with_timing.py`** (script nuevo) es un envoltorio
genérico que se antepone a cualquier comando real:

```
python workflow/scripts/run_with_timing.py \
  --sample-id EC001 --module spades --threads 8 \
  --output results/tables/performance/EC001_spades.tsv \
  -- spades.py -1 R1.fastq.gz -2 R2.fastq.gz -o outdir -t 8 --careful
```

Registra tiempo real, tiempo de CPU, RAM máxima (`resource.getrusage`),
código de salida, fecha e hilos — los seis datos que pide la sección 16 — en
un archivo por (muestra, módulo). Dos detalles importantes, ambos probados:

- **El código de salida del comando real se propaga** como código de salida
  del wrapper: si la herramienta envuelta falla, Snakemake debe seguir
  detectando la regla como fallida, nunca "tragarse" el error solo porque ya
  se registró la métrica.
- **`ru_maxrss` cambia de unidad según el sistema operativo**: kilobytes en
  Linux, bytes en macOS/BSD. Sin ese ajuste, la RAM quedaría sobrestimada
  ~1000× en macOS. Se corrige detectando la plataforma (`platform.system()`).

Se aplicó a las reglas "importantes" (con costo real de cómputo) de todos los
módulos: `download_sample` (que además llenaba un hueco pendiente — nunca se
había creado la regla de Snakemake para la descarga, solo se había probado
por CLI), `fastp`, `spades`, `quast`, `checkm`, `kraken2`, `prokka`,
`amrfinder`. Los pasos livianos (`filter_contigs`, los `parse_*`,
`classify_cephalosporin_genes`, `compare_to_reference`) no se instrumentan,
ya que no son el cuello de botella que esta medición busca identificar.

**`workflow/scripts/combine_performance.py`** (script nuevo) junta todos los
`results/tables/performance/{sample}_{module}.tsv` en tres tablas:
`performance_summary.tsv` (listado largo, muestra+módulo — tal como pide la
sección 16), `performance_by_sample.tsv` (tiempo **total** y RAM **máxima**
por muestra) y `performance_by_module.tsv` (mismo par de métricas por
módulo, para identificar el módulo más costoso — importante: el tiempo se
**suma** entre módulos, pero la RAM se toma como el **máximo**, no la suma,
porque los módulos de una misma muestra no corren simultáneamente).

**`merge_results.py`** se extendió (resolviendo la nota pendiente que quedó
de la parte 17) para incorporar `total_elapsed_seconds` y `peak_max_ram_gb`
por muestra a la tabla maestra — el detalle por módulo se queda en las
tablas de desempeño, mismo principio de "no ensanchar la tabla maestra" ya
aplicado a los genes de AMR.

Probado: el wrapper mide correctamente un comando de ~1s con asignación de
memoria real, y propaga correctamente un código de salida 42 de un comando
que falla a propósito. `combine_performance.py` probado con seis registros
sintéticos de dos muestras: `total_elapsed_seconds` se suma correctamente
entre módulos, `peak_max_ram_gb` toma el máximo (no la suma), y `spades` se
identifica correctamente como el módulo más costoso. La tabla maestra
extendida se probó de punta a punta con el mismo fixture de tres muestras de
la parte 17, confirmando que las nuevas columnas no alteran `final_status`.

### 19. Pruebas de reproducibilidad

**`workflow/scripts/assess_reproducibility.py`** (script nuevo) compara
corridas repetidas de la misma muestra (convención `{base}_run{n}`, ej.
`EC001_run1`, `EC001_run2`, `EC001_run3`) por pares, usando exactamente las
funciones `exact_gene_concordance()` y `jaccard_similarity()` del documento,
más dos comparaciones categóricas adicionales:

- **Genes y alelos**: concordancia exacta + Jaccard sobre el conjunto de
  genes confiables detectados (el símbolo de gen de AMRFinderPlus ya
  codifica el alelo, así que esto cubre ambos requisitos del documento a la vez).
- **Archivo de ensamblaje final**: hash SHA-256 de `contigs.filtered.fasta`
  — idéntico o no entre corridas (SPAdes no es perfectamente determinista
  entre corridas por sus heurísticas multi-hilo, así que esta comparación
  es real, no una formalidad).
- **Estado de clasificación**: coincidencia del `final_status` de la tabla
  maestra entre corridas.

**Límite de arquitectura respetado a propósito:** el documento reserva el
coeficiente de variación **exclusivamente para R** ("R se reservará
únicamente para... Coeficiente de variación"). Por eso este script **no
calcula CV**, ni siquiera como medida secundaria de apoyo — solo usa
comparaciones categóricas/exactas. Los datos numéricos crudos por corrida
(cobertura, tiempo, RAM) ya quedan disponibles en `master_results.tsv` con
un `sample_id` por corrida; R los tomará de ahí (vía la columna `run` de
`validation_input.csv`) para calcular el CV formalmente en la parte 20.

Probado: funciones puras (`exact_gene_concordance`, `jaccard_similarity`,
`parse_replicate_run_id`, `group_replicate_runs`) contra casos conocidos, y
un caso de extremo a extremo con tres corridas sintéticas de `EC001` (dos
idénticas en genes/hash/estado, una tercera con un gen no detectado, hash de
ensamblaje distinto y estado `WARNING` en vez de `PASS`) — las 3
comparaciones por pares resultan correctas: 1 par totalmente concordante,
2 pares discordantes en las tres dimensiones a la vez.

### 20. Estadística en R

**`workflow/scripts/prepare_validation_input.py`** (script nuevo, en Python)
junta `reference_comparison.tsv` con `performance_by_sample.tsv` en
`results/statistics/validation_input.csv`, con exactamente las columnas que
pide el documento: `sample_id, reference_result, pipeline_result, run,
elapsed_seconds, max_ram_gb`. Reutiliza la misma convención de corridas de
reproducibilidad (`{base}_runN`) de la parte 19 para poblar la columna `run`
correctamente; las muestras sin estándar de referencia
(`reference_status="indeterminate"`) se excluyen explícitamente, ya que R
espera un factor binario positive/negative y no hay nada contra qué
compararlas. A partir de aquí, **Python no calcula ninguna estadística** —
solo prepara y limpia los datos.

**`workflow/scripts/run_statistics.R`** es el único script del pipeline que
hace estadística formal. Requiere los paquetes de R `readr`, `caret`, `irr` y
`ggplot2` (`caret` depende además de `stringi`, una dependencia transitiva
que no se resuelve automáticamente en todas las instalaciones). Calcula:

- **Matriz de confusión** (`caret::confusionMatrix`) → `confusion_matrix.txt`
  (texto completo) y `contingency_table.csv` (la tabla 2×2 sola, para
  reutilizar en el reporte HTML sin reparsear texto).
- **Sensibilidad, especificidad y exactitud con intervalos de confianza del
  95%** → `classification_metrics.csv`. `confusionMatrix` ya calcula un IC
  para la exactitud global, pero no para sensibilidad/especificidad por
  separado — se agregó `binom.test()` (Clopper-Pearson, en base R, sin
  dependencias extra) para esas dos.
- **Índice kappa** (`irr::kappa2`) → `kappa.csv`.
- **Coeficiente de variación** de tiempo de ejecución y RAM entre corridas
  repetidas → `cv_execution_time.csv` y `cv_ram_usage.csv` (exactamente la
  función `cv()` del documento, aplicada también a RAM por generalización
  directa, no solo a tiempo).
- **Gráficas**: mapa de calor de la matriz de confusión y barras de CV de
  tiempo por muestra, en `results/statistics/plots/`.

**Límite de arquitectura respetado:** todo el cálculo de CV vive
exclusivamente aquí, en R — Python (partes 18 y 19) solo dejó los datos
crudos por corrida listos, sin calcular ningún CV, tal como exige el
documento ("R se reservará únicamente para... Coeficiente de variación").

Probado de extremo a extremo con datos sintéticos (7 filas: 3 corridas de
reproducibilidad de una muestra + 4 muestras normales, con 1 TP/TN/FP/FN de
cada tipo variado): sensibilidad 0.8 y especificidad 0.5 correctas
aritméticamente (4 TP/1 FN → 0.8; 1 TN/1 FP → 0.5), y — como verificación
cruzada — la kappa calculada independientemente por `irr::kappa2` (0.3)
coincide exactamente con la kappa interna que reporta `caret` en el mismo
análisis. El CV solo se calculó para la muestra con 3 corridas (2.3% en
tiempo, 1.3% en RAM); el resto queda `NA` sin fallar el script (varianza
indefinida con una sola corrida). Corrigiendo un defecto visual encontrado
al revisar el gráfico generado (las celdas con conteo bajo eran casi
invisibles con un gradiente que partía de blanco puro), se cambió la escala
de color para que toda celda sea visible independientemente de su conteo.

### 21. Generación de reportes HTML

**Vacío cerrado antes de empezar:** ningún módulo capturaba la versión de las
herramientas externas (salvo Prokka, ya resuelto en la parte 14), y la
sección 19 exige "Versiones" en cada reporte. Se agregó
**`workflow/scripts/capture_tool_versions.py`** (script nuevo): a diferencia
de las métricas de desempeño (una fila por muestra y módulo), la versión de
una herramienta no cambia entre muestras dentro de una misma corrida, así que
este script se ejecuta **una sola vez por corrida** (no por muestra) y deja
`data/metadata/tool_versions.tsv`. Probado: como ninguna herramienta bio
está instalada en este entorno, todas quedan registradas como
`"not installed"` en vez de hacer fallar el script — la versión es un dato
de trazabilidad, no debe bloquear el pipeline.

**`workflow/templates/sample_report.html.j2`** (carpeta y plantilla nuevas,
no estaban en la estructura fija) es una plantilla Jinja2 con las secciones
exactas que pide la sección 19: identificador/accesiones, calidad de lecturas
(inicial/posterior), cobertura, ensamblaje, taxonomía, completitud/
contaminación, genes detectados (con gráfico embebido), interpretación del
mecanismo, comparación con la referencia, advertencias, desempeño (tiempo y
memoria, total y por módulo) y versiones de herramientas — con un aviso fijo
al inicio de que el reporte no es un diagnóstico clínico.

**`workflow/scripts/generate_report.py`** arma el contexto (uniendo la fila
de la muestra en `master_results.tsv`, sus genes en `amr_classified.tsv`, su
desempeño por módulo, y las versiones de herramientas) y renderiza la
plantilla. Cada reporte es un **HTML autocontenido**: el gráfico de
identidad/cobertura por gen (Matplotlib) se embebe directamente como
imagen en base64, sin archivos `.png` sueltos.

**Regla de alcance respetada literalmente (la más importante de esta
parte):** `build_gene_interpretation_sentences()` nunca genera una
conclusión clínica tipo "Aislado resistente a ceftriaxona". Para cada gen
detectado con confianza, construye una oración con el patrón exacto que pide
el documento: *"Se detectó el determinante blaCTX-M-15, asociado con
beta-lactamasas de espectro extendido (BLEE)..."* — usando la categoría
mecanística ya calculada en la parte 15 (`beta_lactamase_category`), así que
`blaTEM-1` (clasificado como `Other`, no `ESBL`) recibe una interpretación
distinta a `blaCTX-M-15`, en vez de asumir que toda beta-lactamasa detectada
es una BLEE.

Probado de extremo a extremo con `EC001` (3 genes: `blaCTX-M-15` confiable,
`blaTEM-1` confiable, `aac(3)-IId` por debajo del umbral de confianza):

- Verificación automática de que el HTML generado **no contiene** ninguna
  frase de conclusión clínica prohibida (`"aislado resistente"`, etc.).
- `blaCTX-M-15` recibe la interpretación BLEE; `blaTEM-1` recibe la
  interpretación genérica "Other" (no BLEE) — confirma que la regla central
  de la parte 15 se propaga correctamente hasta el reporte final.
- `aac(3)-IId` (bajo el umbral) queda excluido de la interpretación pero sí
  aparece listado en advertencias, sin desaparecer silenciosamente.
- Gráfico embebido en base64 presente en el HTML.
- Probado también el caso sin genes detectados (`EC002`): ambas ramas
  condicionales de la plantilla ("no se detectaron genes...") funcionan.

### 22. Snakefile principal

La validación completa del grafo de dependencias requiere Snakemake
(`pip install snakemake`, versión 9.23.1) — no había sido posible antes,
dado que este archivo es el que finalmente conecta todas las piezas
construidas en las partes anteriores.

**Reglas de agregación finalmente conectadas.** Las reglas `combine_fastp`,
`combine_quast`, `combine_checkm` (+ `checkm_exclusions.tsv`),
`combine_taxonomy` (+ `taxonomy_manual_review.tsv`) y `combine_amr` —
pendientes desde las partes 8 a 15 porque necesitaban `SAMPLES` — se
agregaron a sus archivos `.smk` correspondientes, cada una con
`expand(..., sample=SAMPLES)` sobre la tabla individual de cada muestra.

**Orden correcto de `configfile`/`include`/`SAMPLES`.** El ejemplo del
documento (sección 20) ordena `include:` de las reglas *antes* de definir
`SAMPLES`. Eso no funciona en este pipeline: las reglas `combine_*` nuevas
usan `expand(..., sample=SAMPLES)` directamente en su `input:`, así que
`SAMPLES` debe existir *antes* de que Snakemake evalúe esos archivos. El
`Snakefile` reordena esto: `configfile` → validar muestras y construir
`SAMPLES` → `include:` de todas las reglas → `rule all`.

**Validación de muestras integrada al Snakefile, no como regla aparte.** En
vez de envolver `validate_samples.py` en una regla de Snakemake que produce
un archivo que nadie más usaría, el `Snakefile` importa directamente la
función `validate_samples()` (la misma de la parte 3) y la llama al cargar
el archivo. Si `samples.tsv` tiene un problema, Snakemake falla
inmediatamente al leer el `Snakefile`, antes de intentar construir ningún
grafo con datos no confiables — el script de línea de comandos sigue
disponible para uso manual/CI tal como se probó en la parte 3.

**Dos bugs reales encontrados por el dry-run (`snakemake -n`), no visibles
hasta tener el Snakefile completo:**

1. `combine_performance` declaraba como entrada la *carpeta*
   `results/tables/performance` en vez de la lista de archivos esperados.
   Snakemake no puede saber qué regla "produce" una carpeta suelta —
   necesita depender de archivos concretos. Se corrigió con
   `expand("results/tables/performance/{sample}_{module}.tsv", sample=SAMPLES, module=PERFORMANCE_TRACKED_MODULES)`.
2. **Más importante:** las 7 reglas envueltas con `run_with_timing.py` en la
   parte 18 (`download_sample`, `fastp`, `spades`, `quast`, `checkm`,
   `kraken2`, `prokka`, `amrfinder`) escribían su archivo de desempeño
   dentro del `shell:`, pero **nunca lo declaraban en `output:`**. Snakemake
   solo rastrea archivos declarados explícitamente como salida; un archivo
   que el comando simplemente "escribe de más" es invisible para el grafo de
   dependencias. Se agregó `performance="results/tables/performance/{sample}_{module}.tsv"`
   al `output:` de las 7 reglas, y el `shell:` de cada una ahora referencia
   `{output.performance}` en vez de repetir la ruta a mano (evita que ambas
   copias se desincronicen en el futuro).

Ninguno de los dos bugs habría aparecido probando cada script por separado
como se hizo en las partes anteriores — solo se manifestaron al construir el
grafo completo, que es exactamente para lo que sirve `snakemake -n`.

**`rule all`** pide, por muestra: su reporte HTML; y de forma global: la
tabla maestra, el listado de AMR clasificado, los dos registros de
exclusión/revisión manual (CheckM, Kraken2), y los cinco artefactos de
estadística en R. La anotación (Prokka) queda disponible como regla pero
fuera de `rule all` por defecto: es informativa/complementaria y no forma
parte de los campos exigidos por el reporte individual (sección 19).

**Recomendaciones de `snakemake --lint` no aplicadas, documentadas aquí a
propósito:** el linter sugiere mover los valores de `config[...]` usados
directamente en varios `shell:` hacia la directiva `params:` (mejor
trazabilidad de procedencia), y reemplazar algunos `params` de tipo
"prefijo de ruta" (`outdir`, `output_dir`) por funciones lambda (relevante
para ejecutar en clústeres sin sistema de archivos compartido). Ambas son
recomendaciones de estilo, no errores — confirmado por el dry-run, que
resuelve y ejecuta el grafo completo sin problemas. Se documentan aquí en
vez de aplicarse de forma masiva a las ~15 reglas afectadas, dado que este
proyecto está pensado para correr en una sola máquina (no en un clúster
distribuido) y el beneficio real es marginal frente al costo de tocar tantas
reglas ya probadas.

**Probado:** `snakemake -n --cores 1` (dry-run) contra el `samples.tsv` real
de 3 muestras resuelve el grafo completo: 55 trabajos planeados, cadena
íntegra desde `download_sample` hasta `generate_report` y `run_statistics`
para las 3 muestras. También se verificó con `-p` (imprime los comandos
reales) que los valores de `config.yaml` se sustituyen correctamente en los
comandos generados (ej. `--length_required 50 --qualified_quality_phred 20`
coincide exactamente con `quality.minimum_length` / `quality.minimum_phred`).

### 23. Ambientes Conda

**Simplificación previa necesaria:** `run_with_timing.py` (parte 18) usaba
`pandas`, pero se ejecuta *dentro* del ambiente conda de cada herramienta
que envuelve (fastp, spades, checkm, etc.), no dentro de `python.yaml` —
Snakemake solo activa un ambiente por regla. Eso habría obligado a agregar
`pandas` (y sus dependencias transitivas) a los 8 ambientes de herramientas,
además del propio `python.yaml`. Se reescribió para usar únicamente el
módulo `csv` de la biblioteca estándar; ahora cada ambiente de herramienta
solo necesita agregar `python` (ligero, sin dependencias extra). Re-probado
tras el cambio: mismo comportamiento exacto que antes.

**9 archivos `workflow/envs/*.yaml`** completados (7 existían vacíos + 2
nuevos): `python.yaml` (paquetes de análisis general: pandas, pyyaml,
biopython, jinja2, matplotlib), `fastp.yaml`, `spades.yaml`, `quast.yaml`,
`checkm.yaml`, `kraken2.yaml`, `prokka.yaml`, `amrfinder.yaml` (cada uno con
su herramienta + `python` para `run_with_timing.py`), **`sra_tools.yaml`**
(nuevo — `download_data.py` invoca `fasterq-dump`, que necesita vivir en el
mismo ambiente que `pandas`, no mezclado en `python.yaml`, para respetar el
principio de "un ambiente por herramienta" también para la descarga) y
`r_statistics.yaml` (incluye `r-stringi` explícitamente, la dependencia
transitiva de `caret` que faltó instalar sola en la parte 20).

**Hallazgo importante de compatibilidad de plataforma:** la verificación de
los nombres de paquete contra los índices reales de bioconda (sin instalar,
solo consultando su API) reveló que **QUAST, CheckM (`checkm-genome`) y
Prokka no tienen compilaciones para `osx-arm64`** (Apple Silicon — la
arquitectura de esta Mac y, cada vez más, la más común entre desarrolladores).
Sí existen
para `osx-64` (Intel) y `linux-64`. Esto significa que `snakemake --use-conda`
fallaría al intentar crear esos tres ambientes en esta máquina tal cual.

*Solución práctica estándar para este caso* (bien establecida en la
comunidad de bioinformática para Apple Silicon): forzar a conda a resolver
paquetes `osx-64` y ejecutarlos vía la emulación Rosetta 2, con la variable
de entorno `CONDA_SUBDIR`:

```bash
CONDA_SUBDIR=osx-64 snakemake --use-conda --cores 8 --rerun-incomplete --printshellcmds
```

Esto aplica la emulación a *todos* los ambientes por simplicidad (incluidos
los que sí tienen build nativo arm64, con un costo menor de rendimiento en
esos casos). La alternativa real para producción — y la más alineada con
cómo se plantea usar este pipeline en la práctica (sección 23 del diseño:
"Ejecutar el conjunto independiente de evaluación") — es correr el pipeline
en Linux (contenedor, VM, o un clúster/HPC), donde bioconda tiene
compilaciones completas para las tres herramientas sin necesidad de emulación.

**Probado:** sintaxis YAML válida en los 9 archivos (parseados con
`yaml.safe_load`), y los 8 nombres de paquete bioconda (`fastp`, `spades`,
`quast`, `checkm-genome`, `kraken2`, `prokka`, `ncbi-amrfinderplus`,
`sra-tools`) confirmados como existentes consultando la API de anaconda.org
— no se ejecutó `conda env create` de verdad: los ambientes bioconda más
pesados pueden tardar varios minutos cada uno en resolver, sin necesidad
real de verificarlo en esta etapa del desarrollo.

### 24. Pruebas (unitarias, integración, extremo a extremo, negativas)

La suite completa se construyó en `tests/`, con `pytest` como dependencia
(agregada también a `workflow/envs/python.yaml`). `tests/conftest.py` agrega
`workflow/scripts/` a `sys.path` (los scripts no son un paquete instalado,
son ejecutables independientes) y expone un fixture `repo_root` para que
ningún test dependa del directorio desde el que se invoque `pytest`. Se
agregó `pytest.ini` (`testpaths = tests`) para que `pytest` sin argumentos,
corrido desde la raíz, descubra todo automáticamente.

**106 pruebas, en los 4 niveles que pide la sección 22**, todas reutilizando
los mismos casos que ya se habían validado manualmente a lo largo de las 23
partes anteriores, ahora formalizados con `assert`:

- **Unitarias** (`tests/unit/`, ~96 pruebas): una por script con lógica de
  clasificación — `validate_samples`, `parse_fastp` (incluida la fórmula
  exacta de cobertura del documento), `parse_quast` (incluida la prioridad
  de FAIL sobre WARNING), `parse_checkm`, `parse_kraken2` (incluido el caso
  Shigella), `filter_contigs`, `classify_cephalosporin_genes` (el caso
  crítico TEM/SHV no-BLEE), `compare_to_reference` (incluida la regresión
  del bug de la parte 16), `assess_reproducibility`, `merge_results`.
- **Integración** (`tests/integration/`): verifica el contrato *real* entre
  scripts consecutivos — `parse_amrfinder` → `classify_cephalosporin_genes`
  → `compare_to_reference` encadenados con datos que imitan el formato real
  de AMRFinderPlus, y el ensamblaje completo de la tabla maestra a partir de
  las salidas sintéticas de los 7 módulos.
- **Extremo a extremo** (`tests/e2e/`): una muestra sintética completa
  recorre la cadena entera de scripts en Python — JSON de fastp → report.tsv
  de QUAST → tabla de CheckM → reporte de Kraken2 → tabla de AMRFinderPlus →
  clasificación → comparación de referencia → tabla maestra → reporte HTML
  final — verificando que el HTML resultante contenga el hallazgo genotípico
  correcto y **no** contenga ninguna frase de conclusión clínica prohibida.
- **Negativas** (`tests/negative/`): casos que solo se manifiestan al cruzar
  módulos — una muestra de "especie distinta" o "contaminada" queda
  `EXCLUDED` en la tabla maestra sin desaparecer, pedir el reporte de una
  muestra inexistente falla con un mensaje claro (no un error críptico de
  pandas), y variantes de metadatos incompletos (`samples.tsv` con 0 filas,
  con múltiples columnas faltantes a la vez).

**Alcance declarado explícitamente en los docstrings de las pruebas de
integración y e2e:** ninguna herramienta bioinformática externa (fastp,
SPAdes, QUAST, CheckM, Kraken2, Prokka, AMRFinderPlus) está instalada en
este entorno de desarrollo (y tres de ellas ni siquiera tienen build para
`osx-arm64`, ver parte 23), así que estas pruebas ejercitan la cadena
completa de **código Python propio** con datos que imitan el formato real
de salida de esas herramientas, no las herramientas en sí. Esto es
consistente con cómo se validó cada módulo a lo largo de todo el desarrollo.

**Un hallazgo durante la escritura de las pruebas** (no un bug del pipeline,
sino del primer intento de esta prueba en particular): al construir
`test_low_confidence_detection_excluded...` con una tabla de 2 genes
modificando la cobertura de uno solo, el otro gen (con identidad/cobertura
altas, sin modificar) seguía contando como detección confiable e inflaba el
resultado esperado. Corregido usando una tabla de un solo gen para aislar
exactamente lo que la prueba quería verificar — un recordatorio de que los
fixtures de prueba también pueden tener bugs, no solo el código bajo prueba.

Ejecutar toda la suite: `pytest` (desde la raíz del repositorio).

## Extensiones más allá del diseño original

Tras completar las 24 partes del documento de diseño, se identificaron y
priorizaron algunas extensiones adicionales, no contempladas en el
documento original de 23 secciones:

| # | Parte | Estado |
|---|-------|--------|
| 25 | Segundo motor de AMR (ABricate) y concordancia analítica entre motores | ✅ Hecho |
| 26 | Tipificación de secuencia multilocus (MLST) | ⏳ Pendiente |
| 27 | Descarga de CSV desde el reporte HTML | ⏳ Pendiente |
| 28 | Interfaz web local para análisis ad-hoc (subir FASTQ/FASTA) | ⏳ Pendiente |

### 25. Segundo motor de AMR (ABricate) y concordancia entre motores

Hasta la parte 24, el pipeline comparaba las detecciones de AMRFinderPlus
contra dos referencias distintas: el estándar fenotípico (`compare_to_reference.py`)
y sus propias corridas repetidas (`assess_reproducibility.py`). Ninguna de
las dos mide si una **segunda herramienta de detección independiente**
coincide con AMRFinderPlus sobre el mismo ensamblaje — una señal de
concordancia analítica distinta, y uno de los cuatro pilares que la
descripción general del proyecto declara como objetivo.

**Decisión metodológica central:** comparar genes por **nombre exacto de
alelo** entre AMRFinderPlus (catálogo de referencia de NCBI) y ABricate
(bases de datos CARD + ResFinder) subestimaría artificialmente la
concordancia real, porque cada base de datos nombra los mismos genes de
forma distinta. La comparación se hace a nivel de **familia de gen** (misma
heurística de `derive_gene_family()` ya usada en `parse_amrfinder.py` y
`compare_to_reference.py`, duplicada aquí siguiendo la convención de scripts
autocontenidos del proyecto).

**`workflow/scripts/parse_abricate.py`** (nuevo) normaliza la salida cruda de
ABricate al mismo esquema de columnas que `parse_amrfinder.py`, para que
ambos motores queden en formato directamente comparable. Acepta varios
archivos crudos a la vez (uno por base de datos) y los concatena antes de
normalizar, evitando una regla adicional de Snakemake solo para unir
archivos. Sigue el mismo patrón `parse`/`combine` de archivo-por-muestra.

**`workflow/scripts/compare_amr_engines.py`** (nuevo) calcula, por muestra,
la concordancia (exacta y Jaccard) entre los conjuntos de familias de gen
detectados por cada motor, y prepara además una tabla larga
(`engine_concordance_input.csv`) con el resultado de cada motor como factor
`detected`/`not_detected` por combinación muestra+familia — lista para que R
calcule el índice kappa entre motores. Ningún cálculo estadístico se hace en
Python, mismo principio ya establecido para toda estadística formal del
pipeline.

**`workflow/scripts/compare_engines_statistics.R`** (nuevo, script de R
separado de `run_statistics.R`): calcula el kappa de concordancia
*entre herramientas*, distinto del kappa de `run_statistics.R` (que mide
concordancia contra el estándar de referencia *fenotípico*). Mantenerlos
separados evita mezclar dos preguntas estadísticas distintas en un mismo
script.

**Bug real encontrado y corregido durante las pruebas:** `run_with_timing.py`
(el envoltorio de medición de tiempo/RAM usado por todas las herramientas
externas desde la parte 18) imprimía su resumen de desempeño a `stdout`. Para
la mayoría de las herramientas esto no importaba, porque su resultado real se
escribe en archivos propios (`--output`, `--json`, etc.) y el `stdout`
combinado solo va a un log. Pero ABricate escribe su tabla de resultados
**directamente a stdout** (`abricate ... > resultado.tsv`), así que la regla
`abricate` redirige `stdout` al archivo de resultados — lo que habría
contaminado esa tabla con la línea de estado del envoltorio. Corregido
enviando ese mensaje a `stderr`, comportamiento correcto para una herramienta
de línea de comandos en general (datos por `stdout`, mensajes de estado por
`stderr`), verificado con una prueba dedicada antes y después del cambio.

**Nota de plataforma:** al igual que QUAST/CheckM/Prokka (parte 23),
**`abricate` no tiene compilación para `osx-arm64`** en bioconda (verificado
contra la API de anaconda.org), solo para `osx-64` y `linux-64`. Requiere el
mismo ajuste `CONDA_SUBDIR=osx-64` en esta Mac.

`merge_results.py` incorpora un resumen corto de la concordancia (familias
detectadas por cada motor, concordancia exacta, Jaccard) a la tabla maestra y
al reporte HTML individual, en una sección nueva ("Concordancia entre
motores de AMR"), sin ensanchar la tabla maestra con el detalle completo
(que sigue viviendo en `engine_concordance.tsv`).

Probado: `parse_abricate.py` normaliza correctamente resultados combinados
de CARD + ResFinder, calcula bien el umbral de confianza, y maneja sin
fallar el caso de cero detecciones. `compare_amr_engines.py` distingue
correctamente una discordancia real (un motor detecta, el otro no) de una
concordancia por ausencia (ninguno detecta nada). `compare_engines_statistics.R`
calcula un kappa correcto (0.33) con datos no degenerados, y maneja sin
fallar tanto el caso de varianza indefinida (kappa con muy pocas
observaciones sin variabilidad en un "evaluador", advertencia esperada de
`irr::kappa2`, no un bug) como el caso de tabla completamente vacía. El
grafo completo de Snakemake (`snakemake -n`) resuelve con las 6 ejecuciones
de `abricate` esperadas (2 bases de datos × 3 muestras).

## Estado del roadmap

Las 24 partes del diseño del pipeline están completas. Lo que queda para
llevar esto de "estructura y lógica verificada" a "pipeline ejecutado sobre
datos reales" (pasos 18–22 de la sección 23 del documento original) no es
más código, sino **ejecución real**: instalar los ambientes Conda de verdad
(con el ajuste `CONDA_SUBDIR=osx-64` documentado en la parte 23 si se corre
en esta Mac, o directamente en Linux), descargar un conjunto real de genomas
de *E. coli* desde un repositorio público, correr `snakemake --use-conda`
sobre ellos, revisar los resultados, y solo entonces considerar fijar
versiones de herramientas y publicar. Ningún paso de esos requiere volver a
tocar el código del pipeline salvo que la ejecución real revele un problema
no capturado por las pruebas — que es exactamente el tipo de cosa que estas
pruebas no pueden garantizar al 100%, dado que ninguna corrió contra la
herramienta bioinformática real subyacente.

## Licencia

Este proyecto se distribuye bajo la licencia MIT. El texto completo está en
el archivo [`LICENSE`](LICENSE), en la raíz del repositorio.
