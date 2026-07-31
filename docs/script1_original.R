# -------------------------
# 1️⃣ Librerías
# -------------------------
library(readr)
library(FactoMineR) #Para PCA
library(factoextra) # Para visualización opcional del PCA
library(pheatmap)

# -------------------------
# 2️⃣ Cargar datos
# -------------------------
expr_data <- read_csv("GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")

# Convertir a matriz numérica usando data.matrix de manera segura
expr_matrix <- data.matrix(expr_data[,-1])
rownames(expr_matrix) <- make.names(expr_data[[1]], unique = TRUE)

# -------------------------
# 3️⃣ Heatmap exploratorio (opcional)
# -------------------------
pheatmap(expr_matrix, scale="row", show_rownames=FALSE, show_colnames=TRUE)

# -------------------------
# 4️⃣ PCA con FactoMineR
# -------------------------

# Transponer: filas = muestras, columnas = genes
expr_t <- t(expr_matrix)

# PCA con FactoMineR
# scale.unit = TRUE: escala variables (genes)
pca_res <- PCA(expr_t, scale.unit = TRUE, graph = FALSE)

# -------------------------
# 5️⃣ Ver resultados importantes
# -------------------------

# Varianza explicada por cada componente
pca_res$eig
pca_res$var

# Contribución de cada variable (gen) a los PCs
head(pca_res$var$contrib[,1], 10)  # Top 10 genes que más contribuyen al PC1



# -------------------------
# 6️⃣ Guardar primeros 10 PCs por muestra (para multimodal)
# -------------------------

pc_df <- as.data.frame(pca_res$ind$coord[, 1:5])
pc_df$sample <- rownames(pc_df)
write.csv(pc_df, "GSE96058_FactoMiner_PCA.csv", row.names = FALSE)

# -------------------------
# 7️⃣ Visualización opcional
# -------------------------
fviz_pca_ind(pca_res, repel=TRUE, col.ind="cos2")

# Top 10 genes que más contribuyen a Dim1
top10_genes_dim1 <- head(sort(pca_res$var$contrib[, "Dim.1"], decreasing=TRUE), 10)
top10_genes_dim1
fviz_pca_var(pca_res, select.var=list(contrib=10))
