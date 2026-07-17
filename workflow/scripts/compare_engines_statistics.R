#!/usr/bin/env Rscript
#
# compare_engines_statistics.R
#
# Calcula el indice kappa de concordancia entre los dos motores independientes
# de deteccion de AMR de este pipeline (AMRFinderPlus y ABricate/CARD+ResFinder),
# a partir de la tabla larga que prepara compare_amr_engines.py (Python solo
# prepara y limpia los datos; ningun calculo estadistico se hace ahi).
#
# Este es un analisis de concordancia DISTINTO al de run_statistics.R:
# run_statistics.R mide que tan bien el pipeline reproduce el estandar de
# referencia fenotipico; este script mide que tan bien coinciden dos
# HERRAMIENTAS DE DETECCION independientes entre si, corridas sobre el mismo
# ensamblaje -- una senal de concordancia analitica, no de precision clinica.
#
# Entrada esperada (columnas de engine_concordance_input.csv):
#   sample_id, gene_family, amrfinder_result, abricate_result
# (cada fila es una familia de gen detectada por AL MENOS uno de los dos
# motores en esa muestra; ver compare_amr_engines.py para el detalle)
#
# Uso:
#   Rscript compare_engines_statistics.R [engine_concordance_input.csv] [carpeta_de_salida]

suppressPackageStartupMessages({
  library(readr)
  library(irr)
})

# ------------------------------------------------------------------------
# 1. Argumentos y rutas de salida
# ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_csv_path <- if (length(args) >= 1) args[[1]] else "results/statistics/engine_concordance_input.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "results/statistics"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# 2. Cargar datos y preparar los factores binarios
# ------------------------------------------------------------------------
agreement_data <- read_csv(input_csv_path, show_col_types = FALSE)

if (nrow(agreement_data) == 0) {
  # Ningun motor detecto ninguna familia de gen en ningun lote de muestras:
  # no hay nada que comparar. Se deja constancia en vez de fallar con un
  # error criptico de irr::kappa2 sobre una tabla vacia.
  cat("No hay detecciones de AMR para comparar entre motores (tabla vacia).\n")
  write_csv(
    data.frame(kappa_value = NA, z_statistic = NA, p_value = NA, n_comparisons = 0),
    file.path(output_dir, "engine_concordance_kappa.csv")
  )
  quit(save = "no", status = 0)
}

agreement_data$amrfinder_result <- factor(
  agreement_data$amrfinder_result,
  levels = c("not_detected", "detected")
)
agreement_data$abricate_result <- factor(
  agreement_data$abricate_result,
  levels = c("not_detected", "detected")
)

# ------------------------------------------------------------------------
# 3. Tabla de contingencia (AMRFinderPlus vs. ABricate)
# ------------------------------------------------------------------------
engine_contingency_table <- table(
  amrfinder = agreement_data$amrfinder_result,
  abricate = agreement_data$abricate_result
)
capture.output(
  engine_contingency_table,
  file = file.path(output_dir, "engine_contingency_table.txt")
)

contingency_table_long <- as.data.frame(engine_contingency_table)
write_csv(contingency_table_long, file.path(output_dir, "engine_contingency_table.csv"))

# ------------------------------------------------------------------------
# 4. Indice kappa entre motores
# ------------------------------------------------------------------------
engine_kappa_result <- kappa2(
  data.frame(agreement_data$amrfinder_result, agreement_data$abricate_result)
)
engine_kappa_summary <- data.frame(
  kappa_value = engine_kappa_result$value,
  z_statistic = engine_kappa_result$statistic,
  p_value = engine_kappa_result$p.value,
  n_comparisons = engine_kappa_result$subjects
)
write_csv(engine_kappa_summary, file.path(output_dir, "engine_concordance_kappa.csv"))

cat("Concordancia entre motores de AMR calculada. Resultados en:", output_dir, "\n")
cat("Kappa:", engine_kappa_summary$kappa_value, "\n")
