# scripts/09_c2_differential_expression.R
# ¿QUÉ ES C2? Expresión diferencial del subgrupo topológico C2 (LumA de peor
# pronóstico) frente al resto de LumA (C1 + C4), con limma. Da los genes/biología
# que definen el subgrupo -> saber si es una entidad conocida o algo nuevo.
#
# Requiere: 06_mapper_survival.R (comunidades). Paquete: limma (Bioconductor).
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
if (!requireNamespace("limma", quietly=TRUE)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
  BiocManager::install("limma", update=FALSE, ask=FALSE)
}
suppressMessages(library(limma))

# --- 1. Expresión + comunidades (LumA en C1/C2/C4) ---------------------------
expr <- load_expression(here::here("data","raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
sc <- readr::read_csv(here::here("data","processed","mapper_sample_communities.csv"),
                      show_col_types = FALSE)
sc <- sc[sc$subtype=="LumA" & sc$community %in% c(1,2,4), ]
sc <- sc[sc$sample %in% colnames(expr), ]

E <- expr[, sc$sample]
group <- factor(ifelse(sc$community==2, "C2", "rest"), levels=c("rest","C2"))
cat(sprintf("LumA para DE: %d (C2=%d, resto=%d)\n",
            ncol(E), sum(group=="C2"), sum(group=="rest")))

# --- 2. limma: C2 vs resto ---------------------------------------------------
design <- model.matrix(~ group)                 # coef 2 = C2 vs resto
fit <- eBayes(lmFit(E, design))
tt  <- topTable(fit, coef=2, number=Inf, sort.by="P")
tt$gene <- rownames(tt)
readr::write_csv(tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
                 here::here("results","tables","c2_differential_expression.csv"))

# --- 3. Resumen de significancia ---------------------------------------------
sig <- tt[tt$adj.P.Val < 0.05, ]
cat(sprintf("\nGenes DE (adj.P<0.05): %d  | con |logFC|>0.5: %d  | >1: %d\n",
            nrow(sig), sum(abs(sig$logFC)>0.5), sum(abs(sig$logFC)>1)))

show <- function(df, n=20) print(df[1:min(n,nrow(df)), c("gene","logFC","adj.P.Val")], row.names=FALSE)
up   <- sig[order(-sig$logFC), ]        # más altos en C2
down <- sig[order( sig$logFC), ]        # más bajos en C2
cat("\n=== TOP genes ALTOS en C2 (|logFC| mayor) ===\n"); show(up[abs(up$logFC)>0.5, ], 20)
cat("\n=== TOP genes BAJOS en C2 ===\n");                 show(down[abs(down$logFC)>0.5, ], 20)

# --- 4. Volcano --------------------------------------------------------------
png(here::here("results","figures","c2_volcano.png"), width=1200, height=900, res=150)
plot(tt$logFC, -log10(tt$P.Value), pch=20, cex=0.4, col="grey60",
     xlab="logFC (C2 vs resto LumA)", ylab="-log10(p)", main="DE: C2 vs resto de LumA")
sigp <- tt$adj.P.Val<0.05 & abs(tt$logFC)>0.5
points(tt$logFC[sigp], -log10(tt$P.Value[sigp]), pch=20, cex=0.5,
       col=ifelse(tt$logFC[sigp]>0, "#d73027", "#4575b4"))
abline(v=c(-0.5,0.5), lty=2, col="grey40")
dev.off()

save_session_info(here::here("results","sessionInfo_09_c2_de.txt"))
message("DE de C2 lista: results/tables/c2_differential_expression.csv y figuras.")
