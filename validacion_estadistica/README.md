# Validación estadística

Segunda sección del proyecto, distinta de la herramienta operativa del
pipeline (`Snakefile`, `workflow/`, `webapp/`). Mientras esa primera
sección se enfoca en correr el pipeline sobre una muestra dada y producir
su reporte, esta sección analiza, en conjunto y con muestras reales
aportadas para este fin, qué tan bien concuerdan las predicciones
genotípicas del pipeline con los fenotipos de referencia documentados.

## Estructura

- `muestras/`: metadatos y accesiones de las muestras reales usadas para
  la validación (no se versionan los datos crudos en sí, mismo criterio
  que `data/raw/` en la sección operativa).
- `resultados/`: resultados del pipeline recolectados para cada una de
  esas muestras (tablas, reportes), como entrada al análisis.
- `notebooks/validacion_estadistica.ipynb`: notebook de R (kernel IRkernel)
  con tres apartados — muestras reales, resultados por muestra, y una
  discusión final que integra todos los casos.

## Estado

Estructura creada, aún sin contenido: pendiente de que se aporten las
muestras reales para esta validación.
