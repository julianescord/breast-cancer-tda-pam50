# scripts/15_figures_english.R
# Regenera las 4 figuras del manuscrito con etiquetas/títulos en INGLÉS y correctos
# (Fig 1 = Ball Mapper, no "lente PC1"). Reutiliza datos guardados para que el
# subgrupo C2 sea idéntico al del manuscrito. Mantiene dimensiones originales.
#
# Requiere: 01 (PCs), 03 (etiquetas), 06 (mapper_sample_communities.csv),
#           09 (c2_differential_expression.csv), descarga METABRIC. mappeR, igraph, survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages({library(mappeR); library(igraph); library(survival)})
FIG <- here::here("results","figures")
set.seed(2024)

pcs  <- readr::read_csv(here::here("data","processed","GSE96058_pca_scores.csv"), show_col_types=FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[, c("sample","pam50")]
df   <- dplyr::inner_join(pcs, meta, by="sample")
M    <- as.matrix(df[, grep("^PC", names(df), value=TRUE)])

# ---- Figure 1: Ball Mapper coloured by PAM50 --------------------------------
dmat <- as.matrix(dist(M))
eps  <- as.numeric(quantile(dmat[upper.tri(dmat)], 0.10))
mo <- create_ball_mapper_object(as.data.frame(M), dmat, eps)
verts <- mo[[1]]; edges <- mo[[2]]
npts <- lapply(verts$data, function(s) as.integer(strsplit(s, ",\\s*")[[1]]))
maj  <- sapply(npts, function(i) names(sort(table(factor(df$pam50)[i]), decreasing=TRUE))[1])
pal  <- c(LumA="#2c7fb8", LumB="#7fcdbb", Her2="#f03b20", Basal="#feb24c", Normal="#31a354")
vattr <- data.frame(id=verts$id, size=3+9*sqrt(sapply(npts,length)/max(sapply(npts,length))), color=pal[maj])
g <- graph_from_data_frame(edges[,c("source","target")], directed=FALSE, vertices=vattr)
png(file.path(FIG,"mapper_graph.png"), width=1400, height=1100, res=150)
set.seed(1)
plot(g, layout=layout_with_fr(g), vertex.label=NA, vertex.size=V(g)$size, vertex.color=V(g)$color,
     main="Mapper graph of SCAN-B tumours — colour = majority PAM50 subtype")
legend("topright", legend=names(pal), col=pal, pch=19, bty="n")
dev.off(); message("Fig 1 OK")

# ---- Figure 2: KM within Luminal A by community -----------------------------
sc <- readr::read_csv(here::here("data","processed","mapper_sample_communities.csv"), show_col_types=FALSE)
d2 <- sc[sc$subtype=="LumA" & !is.na(sc$os_days) & sc$os_days>0 & !is.na(sc$os_event), ]
d2 <- d2[d2$community %in% names(which(table(d2$community)>=20)), ]
d2$community <- factor(d2$community)
lr <- survdiff(Surv(os_days, os_event) ~ community, d2); p2 <- 1-pchisq(lr$chisq, length(lr$n)-1)
palc <- rainbow(6)
png(file.path(FIG,"km_within_dominant.png"), width=1100, height=800, res=150)
plot(survfit(Surv(os_days,os_event)~community, d2), col=palc[as.integer(levels(d2$community))], lwd=2,
     xlab="days", ylab="Survival probability",
     main=sprintf("Overall survival within Luminal A by community (log-rank p=%.3g)", p2))
legend("bottomleft", legend=paste0("C", levels(d2$community)),
       col=palc[as.integer(levels(d2$community))], lwd=2, bty="n")
dev.off(); message("Fig 2 OK")

# ---- Figure 3: volcano C2 vs rest of LumA -----------------------------------
de <- readr::read_csv(here::here("results","tables","c2_differential_expression.csv"), show_col_types=FALSE)
sigp <- de$adj.P.Val<0.05 & abs(de$logFC)>0.5
png(file.path(FIG,"c2_volcano.png"), width=1200, height=900, res=150)
plot(de$logFC, -log10(de$P.Value), pch=20, cex=0.4, col="grey60",
     xlab="log2 fold change (C2 vs. rest of Luminal A)", ylab="-log10(p-value)",
     main="Differential expression: C2 vs. remaining Luminal A")
points(de$logFC[sigp], -log10(de$P.Value[sigp]), pch=20, cex=0.5,
       col=ifelse(de$logFC[sigp]>0, "#d73027", "#4575b4"))
abline(v=c(-0.5,0.5), lty=2, col="grey40")
legend("topleft", legend=c("up in C2","down in C2"), col=c("#d73027","#4575b4"), pch=20, bty="n")
dev.off(); message("Fig 3 OK")

# ---- Figure 4: METABRIC KM (C2-like vs rest) --------------------------------
ext <- here::here("data","raw","external")
cl <- data.table::fread(file.path(ext,"metabric_clinical_patient.txt"), skip="PATIENT_ID")[
        , c("PATIENT_ID","CLAUDIN_SUBTYPE","OS_MONTHS","OS_STATUS")]
cl$event <- as.integer(grepl("DECEASED|^1", cl$OS_STATUS)); cl$OS_MONTHS <- suppressWarnings(as.numeric(cl$OS_MONTHS))
ex <- data.table::fread(file.path(ext,"metabric_expression.txt")); sym <- ex[[1]]
Xm <- as.matrix(ex[,-(1:2)]); Xm <- rowsum(Xm, sym, na.rm=TRUE)/as.numeric(table(sym)[rownames(rowsum(Xm,sym))])
Zm <- t(scale(t(Xm)))
mb <- ifelse(de$gene %in% rownames(Xm), de$gene,
             ifelse(gsub(".","-",de$gene,fixed=TRUE) %in% rownames(Xm), gsub(".","-",de$gene,fixed=TRUE), NA))
dd <- de[!is.na(mb),]; dd$mb <- mb[!is.na(mb)]
up <- head(dd$mb[order(-dd$logFC)],100); dn <- head(dd$mb[order(dd$logFC)],100)
c2 <- as.numeric(scale(colMeans(Zm[up,,drop=FALSE],na.rm=TRUE)-colMeans(Zm[dn,,drop=FALSE],na.rm=TRUE)))
sig <- data.frame(PATIENT_ID=colnames(Xm), c2=c2)
m <- merge(cl, sig, by="PATIENT_ID")
m <- m[m$CLAUDIN_SUBTYPE=="LumA" & !is.na(m$OS_MONTHS) & m$OS_MONTHS>0 & !is.na(m$event),]
m$grp <- factor(ifelse(m$c2>median(m$c2),"C2-like","rest"), levels=c("rest","C2-like"))
lr4 <- survdiff(Surv(OS_MONTHS,event)~grp, m); p4 <- 1-pchisq(lr4$chisq,1)
png(file.path(FIG,"metabric_km_c2.png"), width=1100, height=800, res=150)
plot(survfit(Surv(OS_MONTHS,event)~grp, m), col=c("#4575b4","#d73027"), lwd=2,
     xlab="months", ylab="Survival probability",
     main=sprintf("METABRIC · Luminal A: C2-like vs. rest (log-rank p=%.3g)", p4))
legend("bottomleft", legend=levels(m$grp), col=c("#4575b4","#d73027"), lwd=2, bty="n")
dev.off(); message("Fig 4 OK")
message("Las 4 figuras regeneradas en inglés en results/figures/.")
