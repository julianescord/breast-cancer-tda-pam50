# scripts/04_evaluate.R
# Evaluación detallada del SVM ya entrenado: métricas por clase (precisión,
# recall/sensibilidad, especificidad, F1) y matrices de confusión, en TRAIN y TEST.
# Reutiliza el modelo guardado (results/svm_final_workflow.rds); NO re-entrena.
#
# Reconstruye el mismo split de forma determinista (misma semilla, mismo filtro
# de genes) que scripts/02_model_svm.R.
#
# Requiere: haber corrido 02_model_svm.R (genera el .rds) y 03_fetch_labels.R.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
library(tidymodels)
tidymodels_prefer()

# --- 1. Reconstruir los datos exactamente como en 02_model_svm.R -------------
n_top_genes <- as.integer(Sys.getenv("SVM_N_TOP_GENES", "3000"))
expr_path <- here::here(
  "data", "raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
expr <- load_expression(expr_path)
expr <- top_variable_genes(expr, n = n_top_genes)   # filtro determinista
X <- to_model_frame(expr)

meta <- readr::read_csv(here::here("data", "raw", "metadata.csv"),
                        show_col_types = FALSE)
meta <- meta[, c("sample", "pam50", "patient_id")]
meta <- meta[!is.na(meta$pam50) & meta$pam50 != "", ]
dat  <- dplyr::inner_join(meta, X, by = "sample")
dat$outcome <- factor(dat$pam50); dat$pam50 <- NULL

# Mismo split que en 02 (set.seed justo antes de group_initial_split).
set.seed(2024)
split <- group_initial_split(dat, group = patient_id, prop = 0.8)
train <- training(split)
test  <- testing(split)

# --- 2. Cargar el modelo entrenado -------------------------------------------
wf <- readRDS(here::here("results", "svm_final_workflow.rds"))

lev       <- levels(dat$outcome)
prob_cols <- paste0(".pred_", lev)

# --- 3. Evaluar un conjunto: globales + por clase + confusión ----------------
eval_set <- function(df, name) {
  pred <- augment(wf, new_data = df)

  acc <- yardstick::accuracy(pred, outcome, .pred_class)$.estimate
  kap <- yardstick::kap(pred, outcome, .pred_class)$.estimate
  auc <- yardstick::roc_auc(pred, outcome, !!!rlang::syms(prob_cols))$.estimate

  pc  <- per_class_metrics(pred$outcome, pred$.pred_class)
  readr::write_csv(pc, here::here("results", "tables",
                                  paste0("svm_perclass_", name, ".csv")))

  cm <- as.data.frame.matrix(table(Prediction = pred$.pred_class,
                                    Truth = pred$outcome))
  write.csv(cm, here::here("results", "tables",
                           paste0("svm_confusion_", name, ".csv")))

  cat("\n================ ", toupper(name), " ================\n", sep = "")
  cat(sprintf("accuracy=%.3f  balanced_acc=%.3f  macro_F1=%.3f  kappa=%.3f  roc_auc=%.3f\n",
              acc, mean(pc$recall_sensitivity, na.rm = TRUE),
              mean(pc$f1, na.rm = TRUE), kap, auc))
  print(pc, row.names = FALSE)

  data.frame(set = name, accuracy = acc,
             balanced_accuracy = mean(pc$recall_sensitivity, na.rm = TRUE),
             macro_f1 = mean(pc$f1, na.rm = TRUE), kappa = kap, roc_auc = auc)
}

overall <- rbind(eval_set(train, "train"), eval_set(test, "test"))

cat("\n================ RESUMEN TRAIN vs TEST ================\n")
print(overall, row.names = FALSE)
readr::write_csv(overall, here::here("results", "tables",
                                     "svm_overall_train_test.csv"))

save_session_info(here::here("results", "sessionInfo_04_evaluate.txt"))
message("\nEvaluación completa. Tablas por clase y matrices en results/tables/.")
