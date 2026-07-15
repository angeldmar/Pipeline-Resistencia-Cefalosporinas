#!/usr/bin/env Rscript
#
# run_statistics.R
#
# Unico script de este pipeline que hace estadistica formal (sensibilidad,
# especificidad, exactitud, indice kappa, intervalos de confianza,
# coeficiente de variacion y graficas). Toda la preparacion de datos ya se
# hizo en Python (prepare_validation_input.py); este script solo lee
# validation_input.csv y calcula.
#
# Entrada esperada (columnas de validation_input.csv):
#   sample_id, reference_result, pipeline_result, run, elapsed_seconds, max_ram_gb
#
# Uso:
#   Rscript run_statistics.R [validation_input.csv] [carpeta_de_salida]
# (ambos argumentos son opcionales; por defecto usan las rutas estandar del
# pipeline: results/statistics/validation_input.csv y results/statistics)

suppressPackageStartupMessages({
  library(readr)
  library(caret)
  library(irr)
  library(ggplot2)
})

# ------------------------------------------------------------------------
# 1. Argumentos y rutas de salida
# ------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_csv_path <- if (length(args) >= 1) args[[1]] else "results/statistics/validation_input.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "results/statistics"
plots_dir <- file.path(output_dir, "plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# 2. Cargar datos y preparar los factores binarios que espera caret
# ------------------------------------------------------------------------
validation_data <- read_csv(input_csv_path, show_col_types = FALSE)

validation_data$reference_result <- factor(
  validation_data$reference_result,
  levels = c("negative", "positive")
)
validation_data$pipeline_result <- factor(
  validation_data$pipeline_result,
  levels = c("negative", "positive")
)

# ------------------------------------------------------------------------
# 3. Matriz de confusion, sensibilidad, especificidad y exactitud
# ------------------------------------------------------------------------
confusion_result <- confusionMatrix(
  validation_data$pipeline_result,
  validation_data$reference_result,
  positive = "positive"
)

# Salida completa de caret (incluye la tabla, exactitud con su IC del 95%,
# kappa de caret, sensibilidad, especificidad, y demas metricas derivadas).
capture.output(confusion_result, file = file.path(output_dir, "confusion_matrix.txt"))

# Tabla de contingencia (2x2) tambien aparte, en CSV, para que sea facil de
# reutilizar en el reporte HTML sin tener que reparsear el .txt.
contingency_table <- as.data.frame(confusion_result$table)
colnames(contingency_table) <- c("pipeline_result", "reference_result", "count")
write_csv(contingency_table, file.path(output_dir, "contingency_table.csv"))

# Intervalos de confianza (Clopper-Pearson, via binom.test) para
# sensibilidad y especificidad especificamente: confusionMatrix ya calcula
# un IC del 95% para la exactitud global, pero no para sensibilidad y
# especificidad por separado, asi que se calculan aqui con la misma logica.
confusion_counts <- confusion_result$table
true_positive_count <- confusion_counts["positive", "positive"]
false_negative_count <- confusion_counts["negative", "positive"]
true_negative_count <- confusion_counts["negative", "negative"]
false_positive_count <- confusion_counts["positive", "negative"]

compute_rate_with_ci <- function(successes, total) {
  # Devuelve la proporcion y su intervalo de confianza del 95% (Clopper-
  # Pearson). Si no hay observaciones para ese denominador (por ejemplo, cero
  # muestras con referencia positiva), se devuelve NA en vez de fallar.
  if (total == 0) {
    return(list(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  test_result <- binom.test(successes, total)
  list(
    estimate = unname(test_result$estimate),
    ci_lower = test_result$conf.int[1],
    ci_upper = test_result$conf.int[2]
  )
}

sensitivity_stats <- compute_rate_with_ci(true_positive_count, true_positive_count + false_negative_count)
specificity_stats <- compute_rate_with_ci(true_negative_count, true_negative_count + false_positive_count)
accuracy_overall <- confusion_result$overall

classification_metrics <- data.frame(
  metric = c("Sensitivity", "Specificity", "Accuracy"),
  estimate = c(sensitivity_stats$estimate, specificity_stats$estimate, unname(accuracy_overall["Accuracy"])),
  ci_lower = c(sensitivity_stats$ci_lower, specificity_stats$ci_lower, unname(accuracy_overall["AccuracyLower"])),
  ci_upper = c(sensitivity_stats$ci_upper, specificity_stats$ci_upper, unname(accuracy_overall["AccuracyUpper"])),
  n = c(
    true_positive_count + false_negative_count,
    true_negative_count + false_positive_count,
    nrow(validation_data)
  )
)
write_csv(classification_metrics, file.path(output_dir, "classification_metrics.csv"))

# ------------------------------------------------------------------------
# 4. Indice kappa (concordancia entre el estandar de referencia y el
#    resultado del pipeline, corrigiendo por el acuerdo esperado al azar)
# ------------------------------------------------------------------------
kappa_result <- kappa2(
  data.frame(validation_data$reference_result, validation_data$pipeline_result)
)
kappa_summary <- data.frame(
  kappa_value = kappa_result$value,
  z_statistic = kappa_result$statistic,
  p_value = kappa_result$p.value,
  n_samples = kappa_result$subjects
)
write_csv(kappa_summary, file.path(output_dir, "kappa.csv"))

# ------------------------------------------------------------------------
# 5. Coeficiente de variacion (tiempo de ejecucion y RAM) entre corridas
#    repetidas de una misma muestra. El CV es EXCLUSIVO de R en este
#    pipeline; Python solo prepara los datos crudos por corrida.
# ------------------------------------------------------------------------
coefficient_of_variation <- function(x) {
  100 * sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

cv_execution_time <- aggregate(
  elapsed_seconds ~ sample_id,
  data = validation_data,
  FUN = coefficient_of_variation
)
colnames(cv_execution_time) <- c("sample_id", "cv_elapsed_seconds_percent")
write_csv(cv_execution_time, file.path(output_dir, "cv_execution_time.csv"))

cv_ram_usage <- aggregate(
  max_ram_gb ~ sample_id,
  data = validation_data,
  FUN = coefficient_of_variation
)
colnames(cv_ram_usage) <- c("sample_id", "cv_max_ram_gb_percent")
write_csv(cv_ram_usage, file.path(output_dir, "cv_ram_usage.csv"))

# ------------------------------------------------------------------------
# 6. Graficas estadisticas finales
# ------------------------------------------------------------------------

# Mapa de calor de la matriz de confusion.
confusion_matrix_plot <- ggplot(
  contingency_table,
  aes(x = reference_result, y = pipeline_result, fill = count)
) +
  geom_tile(color = "grey40") +
  geom_text(aes(label = count), color = "black", size = 6) +
  # Se evita partir de blanco puro: con pocas muestras, una celda con
  # conteo bajo (ej. 1) quedaria casi invisible sobre el fondo blanco del
  # grafico. "lightyellow" como extremo inferior mantiene cada celda visible
  # independientemente de su conteo.
  scale_fill_gradient(low = "lightyellow", high = "steelblue") +
  labs(
    title = "Matriz de confusion: pipeline vs. estandar de referencia",
    x = "Resultado de referencia",
    y = "Resultado del pipeline",
    fill = "Muestras"
  ) +
  theme_minimal()
ggsave(file.path(plots_dir, "confusion_matrix_heatmap.png"), confusion_matrix_plot, width = 6, height = 5)

# Coeficiente de variacion del tiempo de ejecucion por muestra (solo tiene
# sentido para muestras con mas de una corrida; el resto queda en 0/NA).
cv_time_plot <- ggplot(
  cv_execution_time,
  aes(x = sample_id, y = cv_elapsed_seconds_percent)
) +
  geom_col(fill = "darkorange") +
  labs(
    title = "Coeficiente de variacion del tiempo de ejecucion por muestra",
    x = "Muestra",
    y = "CV (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plots_dir, "cv_execution_time.png"), cv_time_plot, width = 7, height = 5)

cat("Analisis estadistico completo. Resultados en:", output_dir, "\n")
