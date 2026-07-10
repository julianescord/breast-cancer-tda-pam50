# scripts/08_proliferation_cox.R
# VEREDICTO DE NOVEDAD: ¿la comunidad topológica C2 (LumA de peor pronóstico)
# sigue siendo pronóstica tras ajustar por un SCORE DE PROLIFERACIÓN por expresión
# (dato completo en las 3409 muestras, sin los NAs del Ki67-IHC)?
#
# Requiere: 06_mapper_survival.R (comunidades por muestra). Paquete: survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(survival))

# --- 1. Score de proliferación (metagen, firma canónica) ---------------------
expr <- load_expression(here::here("data","raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"))
prolif_genes <- c("ANLN","AURKA","BIRC5","CCNB1","CCNB2","CDC20","CDK1","CDKN3",
                  "CENPF","CEP55","EXO1","FOXM1","KIF2C","MELK","MKI67","NDC80",
                  "NUF2","PTTG1","RRM2","TOP2A","TYMS","UBE2C")
present <- intersect(prolif_genes, rownames(expr))
message(sprintf("Genes de proliferación usados: %d/%d", length(present), length(prolif_genes)))

Z <- t(scale(t(expr[present, ])))            # z-score de cada gen entre muestras
prolif <- colMeans(Z, na.rm = TRUE)          # score por muestra
ps <- data.frame(sample = names(prolif), prolif = as.numeric(scale(prolif)))  # estandarizado (por SD)

# --- 2. Unir con comunidades + grado -----------------------------------------
sc   <- readr::read_csv(here::here("data","processed","mapper_sample_communities.csv"),
                        show_col_types = FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[, c("sample","nhg","ki67_status")]
d <- dplyr::inner_join(sc, meta, by = "sample") |>
     dplyr::inner_join(ps, by = "sample")

d <- d[d$subtype=="LumA" & !is.na(d$os_days) & !is.na(d$os_event) & d$os_days>0, ]
big <- names(which(table(d$community) >= 20))
d <- d[d$community %in% big, ]
d$community <- factor(d$community); d$nhg <- factor(d$nhg)

# --- 3. Validación del score: ¿C2 tiene mayor proliferación? -----------------
cat("\n=== Proliferación media por comunidad (LumA) ===\n")
print(round(tapply(d$prolif, d$community, mean), 3))
cat(sprintf("Correlación score vs Ki67-IHC (donde existe): %.2f\n",
            suppressWarnings(cor(d$prolif, as.numeric(d$ki67_status), use="complete.obs"))))

# --- 4. Cox: comunidad ajustada por grado + proliferación --------------------
dg <- d[!is.na(d$nhg), ]
cat(sprintf("\nLumA analizables (nhg no-NA): %d | proliferación completa\n", nrow(dg)))

cat("\n===== Cox: comunidad + NHG + proliferación =====\n")
m_full <- coxph(Surv(os_days, os_event) ~ community + nhg + prolif, data = dg)
print(round(summary(m_full)$coefficients[, c("exp(coef)","Pr(>|z|)")], 4))

# LRT: ¿la comunidad aporta SOBRE grado + proliferación?
m_base <- coxph(Surv(os_days, os_event) ~ nhg + prolif, data = dg)
lrt <- anova(m_base, m_full)
p_lrt <- lrt$`Pr(>|Chi|)`[2]
cat(sprintf("\nLRT: ¿comunidad aporta sobre NHG + proliferación?  p = %.4g\n", p_lrt))

# --- 5. Veredicto ------------------------------------------------------------
cat("\n===== VEREDICTO DE NOVEDAD =====\n")
if (p_lrt < 0.05) {
  cat(sprintf("La comunidad C2 SIGUE siendo pronóstica tras ajustar por grado Y proliferación\n(LRT p=%.4g): señal INDEPENDIENTE de los factores conocidos -> hallazgo genuino,\nlisto para caracterizar (DE) y validar en cohorte externa.\n", p_lrt))
} else {
  cat(sprintf("La comunidad NO aporta sobre grado + proliferación (LRT p=%.4g): la topología\nre-deriva la proliferación -> validación de biología conocida, no novedad.\n", p_lrt))
}
save_session_info(here::here("results","sessionInfo_08_proliferation_cox.txt"))
