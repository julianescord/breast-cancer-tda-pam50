# scripts/01_eda_pca.R
# Exploración de datos y PCA (SOLO EXPLORATORIO).
#
# ⚠️ IMPORTANTE: el PCA de este script es para explorar/visualizar la estructura
# global de los datos. NO uses estos componentes como entrada del modelo: para
# modelar, el PCA debe ir dentro de la validación cruzada (ver 02_model_svm.R),
# de lo contrario introduces fuga de datos.
#
# Requiere haber ejecutado antes: source("scripts/00_setup.R")
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))

# --- 1. Cargar matriz de expresión (genes x muestras) ------------------------
expr_path <- here::here(
  "data", "raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
expr <- load_expression(expr_path)
message(sprintf("Matriz cargada: %d genes x %d muestras", nrow(expr), ncol(expr)))

# --- 2. Reducir a genes de alta varianza -------------------------------------
# La matriz completa (~30k genes) hace inviable un heatmap y ralentiza el PCA.
# Nos quedamos con los genes más informativos para la exploración.
n_top <- 2000
expr_top <- top_variable_genes(expr, n = n_top)

# --- 3. Heatmap exploratorio (sobre genes top, no la matriz completa) --------
# El clustering de columnas con 3409 muestras es pesado; lo desactivamos por
# defecto. Actívalo (cluster_cols = TRUE) solo si lo necesitas.
png(here::here("results", "figures", "heatmap_top_genes.png"),
    width = 1600, height = 1200, res = 150)
pheatmap::pheatmap(
  expr_top,
  scale         = "row",
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_cols  = FALSE
)
dev.off()

# --- 4. PCA (muestras como observaciones) ------------------------------------
# Transponer: filas = muestras, columnas = genes.
expr_t <- t(expr_top)

# prcomp escala bien y es suficiente para exploración. scale. = TRUE estandariza
# cada gen. (Alternativa FactoMineR::PCA al final del script si la prefieres.)
pca_res <- prcomp(expr_t, center = TRUE, scale. = TRUE)

# Varianza explicada
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
eig_tbl <- data.frame(
  PC              = seq_along(var_explained),
  varianza_pct    = round(100 * var_explained, 3),
  varianza_acum   = round(100 * cumsum(var_explained), 3)
)
head(eig_tbl, 10)
write.csv(eig_tbl,
          here::here("results", "tables", "pca_varianza_explicada.csv"),
          row.names = FALSE)

# --- 5. Scree plot -----------------------------------------------------------
png(here::here("results", "figures", "pca_scree.png"),
    width = 1200, height = 800, res = 150)
barplot(head(eig_tbl$varianza_pct, 20),
        names.arg = head(eig_tbl$PC, 20),
        xlab = "Componente principal", ylab = "% varianza explicada",
        main = "Scree plot (top 20 PCs)")
dev.off()

# --- 6. Guardar scores de las primeras PCs por muestra -----------------------
n_pcs <- 10
pc_df <- as.data.frame(pca_res$x[, seq_len(n_pcs)])
pc_df$sample <- rownames(pc_df)
write.csv(pc_df,
          here::here("data", "processed", "GSE96058_pca_scores.csv"),
          row.names = FALSE)

# --- 7. Genes que más contribuyen a PC1 y PC2 --------------------------------
loadings <- pca_res$rotation
top_pc1 <- head(sort(abs(loadings[, "PC1"]), decreasing = TRUE), 10)
top_pc2 <- head(sort(abs(loadings[, "PC2"]), decreasing = TRUE), 10)
print(top_pc1)
print(top_pc2)

# --- 8. Trazabilidad ---------------------------------------------------------
save_session_info(here::here("results", "sessionInfo_01_eda_pca.txt"))
message("EDA/PCA completado. Salidas en results/ y data/processed/.")

# -----------------------------------------------------------------------------
# Alternativa con FactoMineR (más lenta; úsala sobre expr_t con genes top):
#   library(FactoMineR); library(factoextra)
#   pca_fm <- PCA(expr_t, scale.unit = TRUE, graph = FALSE)
#   fviz_pca_ind(pca_fm, col.ind = "cos2", geom = "point")
#   fviz_pca_var(pca_fm, select.var = list(contrib = 10))
