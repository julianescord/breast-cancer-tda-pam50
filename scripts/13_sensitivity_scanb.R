# scripts/13_sensitivity_scanb.R
# SENSIBILIDAD (SCAN-B): ¿el hallazgo (subgrupo LumA de peor pronóstico, inmuno-frío)
# es robusto a los parámetros del Mapper (eps) y al algoritmo de clustering?
# Para cada configuración: log-rank de supervivencia dentro de LumA entre comunidades,
# y score inmune de la comunidad de PEOR pronóstico (debería ser el más bajo).
#
# Requiere: 01 (PCs), 03 (etiquetas). Paquetes: mappeR, igraph, survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages({library(mappeR); library(igraph); library(survival)})
set.seed(2024)

# --- Datos: PCs + supervivencia + score inmune -------------------------------
pcs  <- readr::read_csv(here::here("data","processed","GSE96058_pca_scores.csv"), show_col_types=FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[, c("sample","pam50","os_days","os_event")]
expr <- load_expression(here::here("data","raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
imm_g <- intersect(c("CD3D","CD3E","CD2","CD8A","CD8B","PTPRC","CD19","MS4A1","CD79A","CD79B",
  "CCL5","CCL19","CXCL9","CXCL10","CXCL13","GZMB","GZMK","IL2RG","IL7R","LCK","CD52","CD48",
  "CD247","CD27"), rownames(expr))
immune <- setNames(colMeans(t(scale(t(expr[imm_g,]))), na.rm=TRUE), colnames(expr))

df <- dplyr::inner_join(pcs, meta, by="sample")
M <- as.matrix(df[, grep("^PC", names(df), value=TRUE)])
dmat <- as.matrix(dist(M))

run_cfg <- function(eps_q, clust) {
  eps <- as.numeric(quantile(dmat[upper.tri(dmat)], eps_q))
  mo <- create_ball_mapper_object(as.data.frame(M), dmat, eps)
  verts <- mo[[1]]; edges <- mo[[2]]
  npts <- lapply(verts$data, function(s) as.integer(strsplit(s, ",\\s*")[[1]]))
  g <- graph_from_data_frame(edges[,c("source","target")], directed=FALSE,
                             vertices=data.frame(id=verts$id))
  cm <- if (clust=="louvain") cluster_louvain(g) else cluster_walktrap(g)
  ncomm <- setNames(as.integer(membership(cm)), names(membership(cm)))[as.character(verts$id)]
  n <- nrow(M); memb <- vector("list", n)
  for (i in seq_along(npts)) for (p in npts[[i]]) memb[[p]] <- c(memb[[p]], ncomm[i])
  sc <- sapply(memb, function(v) if(length(v)) as.integer(names(sort(table(v),decreasing=TRUE))[1]) else NA)

  d <- data.frame(comm=factor(sc), pam=df$pam50, immune=immune[df$sample],
                  os=df$os_days, ev=df$os_event)
  d <- d[d$pam=="LumA" & !is.na(d$os) & d$os>0 & !is.na(d$comm), ]
  big <- names(which(table(d$comm)>=20)); d <- d[d$comm %in% big,]; d$comm <- droplevels(d$comm)
  if (nlevels(d$comm) < 2) return(data.frame(eps_q, clust, ncomm=nlevels(d$comm), logrank_p=NA, worst_immune=NA))
  lr <- survdiff(Surv(os,ev)~comm, d); p <- 1-pchisq(lr$chisq, length(lr$n)-1)
  er <- tapply(d$ev, d$comm, mean); worst <- names(which.max(er))
  data.frame(eps_q, clust, ncomm=nlevels(d$comm), logrank_p=signif(p,3),
             worst_immune=round(mean(d$immune[d$comm==worst]),3))
}

cfgs <- expand.grid(eps_q=c(0.08,0.10,0.12), clust=c("louvain","walktrap"), stringsAsFactors=FALSE)
res <- do.call(rbind, Map(run_cfg, cfgs$eps_q, cfgs$clust))
cat("\n=== SENSIBILIDAD: log-rank dentro de LumA + inmune de la comunidad peor ===\n")
cat("(worst_immune negativo = la comunidad de peor pronóstico es inmuno-fría)\n")
print(res, row.names=FALSE)
readr::write_csv(res, here::here("results","tables","scanb_sensitivity.csv"))
cat(sprintf("\nConfiguraciones con heterogeneidad significativa (p<0.05): %d de %d\n",
            sum(res$logrank_p<0.05, na.rm=TRUE), nrow(res)))
save_session_info(here::here("results","sessionInfo_13_sensitivity.txt"))
