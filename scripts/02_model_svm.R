# scripts/02_model_svm.R
# Clasificador SVM (RBF) del subtipo PAM50 sobre expresión génica, con tidymodels.
# REPRODUCIBLE y SIN FUGA DE DATOS: escalado y PCA se ajustan dentro de cada
# remuestreo (recipe). El único paso fuera de la CV es un filtro de genes por
# varianza NO supervisado (no usa la etiqueta), práctica estándar y necesaria
# por viabilidad computacional (30k genes -> ~3000).
#
# Requiere: source("scripts/00_setup.R") y data/raw/metadata.csv (scripts/03_fetch_labels.R)
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
library(tidymodels)
tidymodels_prefer()
set.seed(2024)

# --- Parámetros (ajusta para pruebas rápidas vs. corrida completa) -----------
# Se pueden sobreescribir con variables de entorno para una prueba rápida, p. ej.:
#   SVM_N_TOP_GENES=300 SVM_CV_FOLDS=3 SVM_GRID_SIZE=3 Rscript scripts/02_model_svm.R
n_top_genes <- as.integer(Sys.getenv("SVM_N_TOP_GENES", "3000"))  # genes por varianza
cv_folds    <- as.integer(Sys.getenv("SVM_CV_FOLDS",    "5"))     # folds de CV
grid_size   <- as.integer(Sys.getenv("SVM_GRID_SIZE",   "15"))    # combinaciones de hiperparámetros

# --- 1. Expresión (genes x muestras) -> pre-filtro varianza -> muestras x genes
expr_path <- here::here(
  "data", "raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
expr <- load_expression(expr_path)
expr <- top_variable_genes(expr, n = n_top_genes)  # filtro NO supervisado
X <- to_model_frame(expr)                          # muestras x genes + `sample`

# --- 2. Etiquetas (subtipo PAM50) --------------------------------------------
labels_path <- here::here("data", "raw", "metadata.csv")
if (!file.exists(labels_path)) {
  stop("Falta data/raw/metadata.csv. Ejecuta antes scripts/03_fetch_labels.R")
}
meta <- readr::read_csv(labels_path, show_col_types = FALSE)

# Solo columnas necesarias para NO contaminar los predictores con datos clínicos:
# la respuesta (pam50) y el identificador de agrupación (patient_id).
meta <- meta[, c("sample", "pam50", "patient_id")]
meta <- meta[!is.na(meta$pam50) & meta$pam50 != "", ]

dat <- dplyr::inner_join(meta, X, by = "sample")
dat$outcome <- factor(dat$pam50)
dat$pam50 <- NULL
message(sprintf("Muestras para modelar: %d | clases: %s",
                nrow(dat), paste(levels(dat$outcome), collapse = ", ")))
print(table(dat$outcome))

id_cols <- c("sample", "patient_id")

# --- 3. Partición train/test agrupada por paciente (réplicas juntas) ---------
# Agrupar por patient_id evita que F#### y F####repl caigan en lados distintos.
split <- group_initial_split(dat, group = patient_id, prop = 0.8)
train <- training(split)
test  <- testing(split)

# --- 4. Remuestreo para tuning (CV agrupada) ---------------------------------
folds <- group_vfold_cv(train, group = patient_id, v = cv_folds)

# --- 5. Receta: preprocesado SIN fuga (se reajusta en cada fold) -------------
rec <- recipe(outcome ~ ., data = train) |>
  update_role(dplyr::all_of(id_cols), new_role = "id") |>
  step_zv(all_numeric_predictors()) |>
  step_normalize(all_numeric_predictors()) |>
  step_pca(all_numeric_predictors(), num_comp = tune())

# --- 6. Modelo SVM (RBF) con hiperparámetros a tunear ------------------------
svm_spec <- svm_rbf(cost = tune(), rbf_sigma = tune()) |>
  set_engine("kernlab") |>
  set_mode("classification")

wf <- workflow() |>
  add_recipe(rec) |>
  add_model(svm_spec)

# --- 7. Malla de hiperparámetros ---------------------------------------------
params <- extract_parameter_set_dials(wf) |>
  update(
    num_comp  = num_comp(range = c(10L, 50L)),
    cost      = cost(range = c(-2, 6)),      # log2: evita cost tan bajo que colapse a la clase mayoritaria
    rbf_sigma = rbf_sigma(range = c(-4, 0))  # log10
  )
grid <- grid_space_filling(params, size = grid_size, type = "latin_hypercube")

# --- 8. Métricas multiclase (clases desbalanceadas) --------------------------
# roc_auc: Hand-Till. bal_accuracy/kap robustas al desbalanceo. sens/spec/
# precision/recall/f_meas se promedian (macro) en multiclase. El desglose POR
# CLASE se obtiene en scripts/04_evaluate.R.
cls_metrics <- metric_set(roc_auc, accuracy, bal_accuracy, kap,
                          sens, spec, precision, recall, f_meas)

# --- 9. Tuning por validación cruzada ----------------------------------------
tuned <- tune_grid(
  wf, resamples = folds, grid = grid, metrics = cls_metrics,
  control = control_grid(save_pred = TRUE, verbose = TRUE)
)

autoplot(tuned)
ggplot2::ggsave(here::here("results", "figures", "svm_tuning.png"),
                width = 9, height = 6, dpi = 150)

# "Evolución" del ajuste = rendimiento (CV) de cada combinación de hiperparámetros.
# (El SVM no tiene epochs; esto es el análogo: la búsqueda por validación cruzada.)
readr::write_csv(collect_metrics(tuned),
                 here::here("results", "tables", "svm_tuning_metrics.csv"))

best <- select_best(tuned, metric = "roc_auc")
print(best)
# Rendimiento en CV del mejor modelo (estimación honesta durante el entrenamiento).
cv_best <- collect_metrics(tuned) |>
  dplyr::semi_join(best, by = c("cost", "rbf_sigma", "num_comp"))
readr::write_csv(cv_best,
                 here::here("results", "tables", "svm_cv_best_metrics.csv"))

# --- 10. Ajuste final y evaluación en el test (una sola vez) ------------------
final_wf  <- finalize_workflow(wf, best)
final_fit <- last_fit(final_wf, split, metrics = cls_metrics)

test_metrics <- collect_metrics(final_fit)
print(test_metrics)
readr::write_csv(test_metrics,
                 here::here("results", "tables", "svm_test_metrics.csv"))

preds <- collect_predictions(final_fit)
cm <- conf_mat(preds, truth = outcome, estimate = .pred_class)
print(cm)
png(here::here("results", "figures", "svm_confusion_matrix.png"),
    width = 1000, height = 900, res = 150)
autoplot(cm, type = "heatmap")
dev.off()

# --- 11. Guardar modelo final y trazabilidad ---------------------------------
saveRDS(extract_workflow(final_fit),
        here::here("results", "svm_final_workflow.rds"))
save_session_info(here::here("results", "sessionInfo_02_model_svm.txt"))
message("Modelado SVM completado. Métricas y modelo en results/.")

# -----------------------------------------------------------------------------
# Notas (ver CLAUDE.md):
# - Para estimaciones aún más honestas: CV anidada (nested_cv).
# - Comparar kernels: svm_linear() / svm_poly().
# - Reportar según DOME (y TRIPOD+AI si se usa un desenlace clínico).
