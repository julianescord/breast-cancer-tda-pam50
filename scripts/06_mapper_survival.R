# scripts/06_mapper_survival.R
# ¿La estructura topológica (Mapper) estratifica la SUPERVIVENCIA dentro de un
# subtipo PAM50, más allá del propio subtipo?  (Patrón del descubrimiento c-MYB+.)
#
# Rigor anti-circular: los subgrupos se definen por TOPOLOGÍA (comunidades del
# grafo Mapper, sin usar supervivencia); la supervivencia se compara DESPUÉS
# (Kaplan-Meier + log-rank).
#
# Requiere: 01_eda_pca.R, 03_fetch_labels.R. Paquetes: mappeR, igraph, survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages({library(mappeR); library(igraph)})
if (!requireNamespace("survival", quietly=TRUE)) install.packages("survival", type="source")
suppressMessages(library(survival))
set.seed(2024)

eps_q <- as.numeric(Sys.getenv("MAP_EPS_Q", "0.10"))

# --- 1. Datos: PCs + subtipo + supervivencia ---------------------------------
pcs  <- readr::read_csv(here::here("data","processed","GSE96058_pca_scores.csv"),
                        show_col_types = FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[
          , c("sample","pam50","os_days","os_event")]
df <- dplyr::inner_join(pcs, meta, by = "sample")
lab      <- factor(df$pam50)
os_days  <- suppressWarnings(as.numeric(df$os_days))
os_event <- suppressWarnings(as.integer(as.character(df$os_event)))
M <- as.matrix(df[, grep("^PC", names(df), value = TRUE)])
message(sprintf("Muestras: %d | con supervivencia: %d",
                nrow(M), sum(!is.na(os_days) & !is.na(os_event))))

# --- 2. Ball Mapper ----------------------------------------------------------
dmat <- as.matrix(dist(M))
eps  <- as.numeric(quantile(dmat[upper.tri(dmat)], eps_q))
mo <- create_ball_mapper_object(as.data.frame(M), dmat, eps)
verts <- mo[[1]]; edges <- mo[[2]]
node_pts <- lapply(verts$data, function(s) as.integer(strsplit(s, ",\\s*")[[1]]))
message(sprintf("Mapper: %d nodos, %d aristas (eps=%.2f)", nrow(verts), nrow(edges), eps))

# --- 3. Comunidades topológicas (SIN mirar supervivencia) --------------------
g <- graph_from_data_frame(edges[, c("source","target")], directed = FALSE,
                           vertices = data.frame(id = verts$id))
comm <- cluster_louvain(g)
node_comm <- setNames(as.integer(membership(comm)), names(membership(comm)))
node_comm_vec <- node_comm[as.character(verts$id)]           # comunidad por nodo
message(sprintf("Comunidades topológicas detectadas: %d", max(node_comm_vec, na.rm=TRUE)))

# Asignar cada muestra a una comunidad (mayoría de los nodos que la contienen)
n <- nrow(M); samp_comm <- integer(n)
memb <- vector("list", n)
for (i in seq_along(node_pts)) for (p in node_pts[[i]]) memb[[p]] <- c(memb[[p]], node_comm_vec[i])
samp_comm <- sapply(memb, function(v) if (length(v)) as.integer(names(sort(table(v), decreasing=TRUE))[1]) else NA)

# Guardar asignación (para análisis posteriores, p. ej. Cox). ANTI-PSEUDO-REPLICACIÓN:
# el Mapper se construyó con las 3409 muestras (incl. réplicas técnicas), pero para la
# INFERENCIA de supervivencia dejamos UNA fila por paciente (F####repl es réplica de
# F####; mismo tumor y mismo desenlace no pueden contar dos veces).
assign_df <- data.frame(sample = df$sample, subtype = as.character(lab),
                        community = samp_comm, os_days = os_days, os_event = os_event,
                        stringsAsFactors = FALSE)
assign_df$patient_id <- sub("repl$", "", assign_df$sample)
sc_d <- drop_technical_replicates(assign_df, "sample")   # 1 fila / paciente
message(sprintf("Filas: %d muestras -> %d pacientes (quitadas %d réplicas técnicas)",
                nrow(assign_df), nrow(sc_d), nrow(assign_df) - nrow(sc_d)))
readr::write_csv(sc_d, here::here("data","processed","mapper_sample_communities.csv"))

# --- 4. Resumen por comunidad (composición de subtipo + supervivencia) -------
comm_summary <- do.call(rbind, lapply(sort(unique(na.omit(sc_d$community))), function(cc) {
  in_c <- which(sc_d$community == cc)
  sv   <- !is.na(sc_d$os_days[in_c]) & !is.na(sc_d$os_event[in_c])
  top  <- sort(table(sc_d$subtype[in_c]), decreasing = TRUE)
  data.frame(comunidad = cc, n = length(in_c),
             subtipo_dom = names(top)[1],
             pct_dom = round(100*top[1]/length(in_c)),
             n_surv = sum(sv),
             tasa_eventos = round(mean(sc_d$os_event[in_c][sv]), 3),
             mediana_dias = round(median(sc_d$os_days[in_c][sv])))
}))
readr::write_csv(comm_summary, here::here("results","tables","mapper_community_summary.csv"))
cat("\n=== Comunidades topológicas (1 fila/paciente) ===\n"); print(comm_summary, row.names = FALSE)

# --- 5. Supervivencia global entre comunidades (log-rank) --------------------
sd <- sc_d[, c("os_days","os_event","community","subtype")]
sd$community <- factor(sd$community)
sd <- sd[!is.na(sd$os_days) & !is.na(sd$os_event) & !is.na(sd$community) & sd$os_days > 0, ]
lr_global <- survdiff(Surv(os_days, os_event) ~ community, data = sd)
p_global <- 1 - pchisq(lr_global$chisq, length(lr_global$n) - 1)
cat(sprintf("\nLog-rank GLOBAL entre comunidades: p = %.3g\n", p_global))

pal <- rainbow(max(node_comm_vec, na.rm=TRUE))
png(here::here("results","figures","km_communities.png"), width=1100, height=800, res=150)
plot(survfit(Surv(os_days, os_event) ~ community, data = sd), col = pal, lwd = 2,
     xlab = "días", ylab = "Supervivencia", main = sprintf("KM por comunidad topológica (log-rank p=%.3g)", p_global))
legend("bottomleft", legend = paste("C", levels(sd$community)), col = pal, lwd = 2, bty = "n")
dev.off()

# --- 6. ¿Estratifica DENTRO del subtipo dominante? (el hallazgo clave) --------
dom <- names(sort(table(sc_d$subtype), decreasing = TRUE))[1]   # LumA
sub <- sd[sd$subtype == dom, ]
sub$community <- droplevels(sub$community)
keep <- names(which(table(sub$community) >= 20))       # comunidades con muestras suficientes
sub <- sub[sub$community %in% keep, ]; sub$community <- droplevels(sub$community)
cat(sprintf("\n=== Dentro de %s: %d muestras, %d comunidades (n>=20) ===\n",
            dom, nrow(sub), nlevels(sub$community)))
if (nlevels(sub$community) >= 2) {
  lr_sub <- survdiff(Surv(os_days, os_event) ~ community, data = sub)
  p_sub <- 1 - pchisq(lr_sub$chisq, length(lr_sub$n) - 1)
  cat(sprintf("Log-rank DENTRO de %s entre comunidades: p = %.3g\n", dom, p_sub))
  png(here::here("results","figures","km_within_dominant.png"), width=1100, height=800, res=150)
  plot(survfit(Surv(os_days, os_event) ~ community, data = sub),
       col = pal[as.integer(levels(sub$community))], lwd = 2,
       xlab = "días", ylab = "Supervivencia",
       main = sprintf("KM dentro de %s por comunidad (log-rank p=%.3g)", dom, p_sub))
  legend("bottomleft", legend = paste("C", levels(sub$community)),
         col = pal[as.integer(levels(sub$community))], lwd = 2, bty = "n")
  dev.off()
} else cat("No hay >=2 comunidades con suficientes muestras dentro del subtipo dominante.\n")

# --- 7. Grafo coloreado por comunidad ----------------------------------------
V(g)$size  <- 3 + 9*sqrt(sapply(node_pts, length)/max(sapply(node_pts, length)))
V(g)$color <- pal[node_comm_vec]
png(here::here("results","figures","mapper_communities.png"), width=1400, height=1100, res=150)
set.seed(1)
plot(g, layout = layout_with_fr(g), vertex.label = NA,
     main = "Comunidades topológicas del Mapper (color = comunidad)")
legend("topright", legend = paste("C", seq_len(max(node_comm_vec))), col = pal, pch = 19, bty = "n")
dev.off()

save_session_info(here::here("results","sessionInfo_06_mapper_survival.txt"))
message("Análisis de supervivencia listo: KM y grafos en results/figures/, tabla en results/tables/.")
