# R/functions.R
# Funciones reutilizables del proyecto. Sin efectos secundarios (no leen/escriben
# fuera de lo que reciben por argumento, no cargan librerías globales).
# Cargar con: source(here::here("R", "functions.R"))

#' Carga la matriz de expresión de GSE96058
#'
#' @param path Ruta al .csv.gz (filas = genes, columnas = muestras).
#' @return Matriz numérica con genes en filas y muestras en columnas.
load_expression <- function(path) {
  stopifnot(file.exists(path))
  # fread lee gzip directamente y es rápido con archivos grandes.
  dt <- data.table::fread(path)
  gene_ids <- dt[[1]]
  m <- as.matrix(dt[, -1L])
  rownames(m) <- make.names(gene_ids, unique = TRUE)
  storage.mode(m) <- "double"
  m
}

#' Selecciona los N genes de mayor varianza
#'
#' Reduce la dimensionalidad antes de PCA/heatmap para que sea tratable.
#' OJO: si se usa como filtro para modelar, debe hacerse DENTRO de la CV
#' (ver step_* en tidymodels), no sobre la matriz completa.
#'
#' @param expr Matriz genes x muestras.
#' @param n Número de genes a conservar.
#' @return Matriz genes x muestras con los n genes más variables.
top_variable_genes <- function(expr, n = 2000) {
  vars <- matrixStats::rowVars(expr, na.rm = TRUE)
  keep <- order(vars, decreasing = TRUE)[seq_len(min(n, nrow(expr)))]
  expr[keep, , drop = FALSE]
}

#' Transpone a formato de modelado: muestras en filas, genes en columnas
#'
#' @param expr Matriz genes x muestras.
#' @return data.frame muestras x genes con columna `sample`.
to_model_frame <- function(expr) {
  df <- as.data.frame(t(expr))
  df$sample <- rownames(df)
  df
}

#' Guarda sessionInfo() para trazabilidad reproducible
#'
#' @param path Ruta del archivo de texto de salida.
save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
  invisible(path)
}

#' Métricas por clase (one-vs-rest) a partir de la verdad y la predicción
#'
#' Calcula soporte, precisión, recall (= sensibilidad), especificidad y F1 para
#' cada clase, derivadas de la matriz de confusión. Sirve para multiclase.
#'
#' @param truth Factor con las clases verdaderas.
#' @param pred_class Vector/factor con las clases predichas.
#' @return data.frame con una fila por clase.
per_class_metrics <- function(truth, pred_class) {
  lev <- levels(truth)
  cm  <- table(Prediction = factor(pred_class, levels = lev),
               Truth      = factor(truth,      levels = lev))
  rows <- lapply(lev, function(k) {
    TP <- cm[k, k]
    FP <- sum(cm[k, ]) - TP
    FN <- sum(cm[, k]) - TP
    TN <- sum(cm) - TP - FP - FN
    prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_   # sensibilidad / recall
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_   # especificidad
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0)
              2 * prec * rec / (prec + rec) else NA_real_
    data.frame(class = k, support = sum(cm[, k]),
               precision = round(prec, 4), recall_sensitivity = round(rec, 4),
               specificity = round(spec, 4), f1 = round(f1, 4))
  })
  do.call(rbind, rows)
}
