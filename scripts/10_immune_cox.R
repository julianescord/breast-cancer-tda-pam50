# scripts/10_immune_cox.R
# ¿C2 va MÁS ALLÁ del eje inmune? Cox dentro de LumA ajustando por grado,
# proliferación y un SCORE DE INFILTRACIÓN INMUNE (expresión, dato completo).
#   - Si C2 sobrevive al ajuste inmune -> señal más allá de la infiltración inmune.
#   - Si desaparece -> C2 = luminal inmuno-frío (eje conocido, hallado por topología).
#
# Requiere: 06_mapper_survival.R. Paquete: survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(survival))

expr <- load_expression(here::here("data","raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))

metagene <- function(genes) {
  g <- intersect(genes, rownames(expr))
  as.numeric(scale(colMeans(t(scale(t(expr[g, ]))), na.rm = TRUE)))
}
prolif_g <- c("ANLN","AURKA","BIRC5","CCNB1","CCNB2","CDC20","CDK1","CDKN3","CENPF",
              "CEP55","EXO1","FOXM1","KIF2C","MELK","MKI67","NDC80","NUF2","PTTG1",
              "RRM2","TOP2A","TYMS","UBE2C")
immune_g <- c("CD3D","CD3E","CD2","CD8A","CD8B","PTPRC","CD19","MS4A1","CD79A","CD79B",
              "CCL5","CCL19","CXCL9","CXCL10","CXCL13","GZMB","GZMK","IL2RG","IL7R",
              "LCK","CD52","CD48","CD247","CD27")
scores <- data.frame(sample = colnames(expr),
                     prolif = metagene(prolif_g),
                     immune = metagene(immune_g))

sc   <- readr::read_csv(here::here("data","processed","mapper_sample_communities.csv"),
                        show_col_types = FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[, c("sample","nhg")]
d <- dplyr::inner_join(sc, meta, by="sample") |> dplyr::inner_join(scores, by="sample")
d <- d[d$subtype=="LumA" & !is.na(d$os_days) & !is.na(d$os_event) & d$os_days>0, ]
d <- d[d$community %in% names(which(table(d$community) >= 20)), ]
d$community <- factor(d$community); d$nhg <- factor(d$nhg)

cat("=== Score inmune medio por comunidad (LumA) ===\n")
print(round(tapply(d$immune, d$community, mean), 3))
cat("(C2 debería ser el más bajo = inmuno-frío)\n")

# (a) comunidad + inmune
cat("\n===== (a) Cox: comunidad + inmune =====\n")
print(round(summary(coxph(Surv(os_days,os_event) ~ community + immune, d))$coefficients[,c("exp(coef)","Pr(>|z|)")],4))

# (b) modelo completo: grado + proliferación + inmune
dg <- d[!is.na(d$nhg), ]
cat(sprintf("\n===== (b) Cox: comunidad + NHG + proliferación + inmune (n=%d) =====\n", nrow(dg)))
m_full <- coxph(Surv(os_days,os_event) ~ community + nhg + prolif + immune, dg)
print(round(summary(m_full)$coefficients[,c("exp(coef)","Pr(>|z|)")],4))

m_base <- coxph(Surv(os_days,os_event) ~ nhg + prolif + immune, dg)
p_lrt <- anova(m_base, m_full)$`Pr(>|Chi|)`[2]
cat(sprintf("\nLRT: ¿comunidad aporta sobre grado+proliferación+inmune?  p = %.4g\n", p_lrt))

cat("\n===== VEREDICTO =====\n")
if (p_lrt < 0.05) cat(sprintf("C2 SIGUE pronóstica tras ajustar tambien por inmune (p=%.4g): va MÁS ALLÁ del eje inmune.\n", p_lrt)) else cat(sprintf("C2 se explica por el eje inmune (p=%.4g): la topología halló, no supervisadamente, el subgrupo luminal inmuno-frío.\n", p_lrt))
save_session_info(here::here("results","sessionInfo_10_immune_cox.txt"))
