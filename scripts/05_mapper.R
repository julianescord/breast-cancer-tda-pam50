# scripts/05_mapper.R
# TDA / Mapper: grafo de la estructura de las 3409 muestras para ver si recupera
# los subtipos PAM50 y si aparecen sub-ramas (posible biología nueva, tipo c-MYB+).
#
# Espacio métrico: 10 PCs (data/processed/GSE96058_pca_scores.csv).
# Lente/filtro 1D: PC1 (eje ER/luminal-basal). Nodos por subtipo PAM50 mayoritario.
#
# Requiere: 01_eda_pca.R (PCs), 03_fetch_labels.R (etiquetas). Paquetes: mappeR, igraph.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages({library(mappeR); library(igraph)})
set.seed(2024)

n_bins      <- as.integer(Sys.getenv("MAP_BINS",    "10"))  # intervalos en la lente
pct_overlap <- as.numeric(Sys.getenv("MAP_OVERLAP", "50"))  # solape %

# --- 1. Datos: PCs + etiquetas -----------------------------------------------
pcs  <- readr::read_csv(here::here("data","processed","GSE96058_pca_scores.csv"),
                        show_col_types = FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[, c("sample","pam50")]
df   <- dplyr::inner_join(pcs, meta, by = "sample")
lab  <- factor(df$pam50)

pc_cols <- grep("^PC", names(df), value = TRUE)
M <- as.matrix(df[, pc_cols])
message(sprintf("Muestras: %d | PCs: %s", nrow(M), paste(pc_cols, collapse=",")))

# --- 2. Construir el objeto Mapper (1D lente PC1  ó  Ball Mapper) -------------
dmat   <- as.matrix(dist(M))
method <- Sys.getenv("MAP_METHOD", "1d")     # "1d" | "ball"
if (method == "ball") {
  epsq <- as.numeric(Sys.getenv("MAP_EPS_Q", "0.12"))
  eps  <- as.numeric(quantile(dmat[upper.tri(dmat)], epsq))
  mo <- create_ball_mapper_object(as.data.frame(M), dmat, eps)
  message(sprintf("Ball Mapper: eps=%.2f (cuantil %.2f de distancias)", eps, epsq))
} else {
  cover <- create_width_balanced_cover(min(M[,1]), max(M[,1]), n_bins, pct_overlap)
  mo <- create_1D_mapper_object(as.data.frame(M), dmat, M[,1], cover,
          global_hierarchical_clusterer("single", dmat))
  message(sprintf("Mapper 1D: lente PC1, %d bins, %.0f%% solape", n_bins, pct_overlap))
}
verts <- mo[[1]]; edges <- mo[[2]]
message(sprintf("Grafo Mapper: %d nodos, %d aristas | cols: %s",
                nrow(verts), nrow(edges), paste(names(verts), collapse=",")))

# --- 3. Subtipo mayoritario y pureza por nodo --------------------------------
node_pts <- lapply(verts$data, function(s) as.integer(strsplit(s, ",\\s*")[[1]]))
summ <- do.call(rbind, lapply(seq_along(node_pts), function(i) {
  tb <- sort(table(lab[node_pts[[i]]]), decreasing = TRUE)
  data.frame(id = verts$id[i], size = length(node_pts[[i]]),
             majority = names(tb)[1], purity = round(tb[1]/length(node_pts[[i]]), 3))
}))
readr::write_csv(summ, here::here("results","tables","mapper_nodes.csv"))

cat(sprintf("\nPureza media de nodos: %.3f  (1 = nodo de un solo subtipo)\n", mean(summ$purity)))
cat(sprintf("Nodos 'puros' (pureza>0.7): %d de %d\n", sum(summ$purity>0.7), nrow(summ)))
cat("Subtipo mayoritario por nodo:\n"); print(table(summ$majority))

# --- 4. Grafo coloreado + guardado -------------------------------------------
pal <- c(LumA="#2c7fb8", LumB="#7fcdbb", Her2="#f03b20", Basal="#feb24c", Normal="#31a354")
vattr <- data.frame(id = verts$id,
                    size = 3 + 9*sqrt(summ$size/max(summ$size)),
                    color = pal[summ$majority])
g <- graph_from_data_frame(edges[, c("source","target")], directed = FALSE, vertices = vattr)

png(here::here("results","figures","mapper_graph.png"), width=1400, height=1100, res=150)
set.seed(1)
plot(g, layout = layout_with_fr(g), vertex.label = NA,
     vertex.size = V(g)$size, vertex.color = V(g)$color,
     main = "Mapper GSE96058 (lente PC1) — color = subtipo PAM50 mayoritario")
legend("topright", legend = names(pal), col = pal, pch = 19, bty = "n")
dev.off()

save_session_info(here::here("results","sessionInfo_05_mapper.txt"))
message("Mapper listo: results/figures/mapper_graph.png y results/tables/mapper_nodes.csv")
