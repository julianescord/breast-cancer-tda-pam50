# scripts/07_cox_luma.R
# TEST DECISIVO DE NOVEDAD: ¿la comunidad topológica sigue siendo pronóstica
# DENTRO de LumA tras ajustar por factores clínicos conocidos (grado NHG, Ki67)?
#   - Si la comunidad se mantiene significativa -> señal pronóstica INDEPENDIENTE (novedad).
#   - Si no -> la topología re-deriva proliferación/grado (validación, no descubrimiento).
#
# Requiere: 06_mapper_survival.R (genera data/processed/mapper_sample_communities.csv).
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(survival))

# --- 1. Comunidades por muestra + covariables clínicas -----------------------
sc <- readr::read_csv(here::here("data","processed","mapper_sample_communities.csv"),
                      show_col_types = FALSE)
meta <- data.table::fread(here::here("data","raw","metadata.csv"))[
          , c("sample","nhg","ki67_status","er_status","her2_status")]
d <- dplyr::inner_join(sc, meta, by = "sample")

# LumA, con supervivencia válida y en las comunidades grandes (n>=20)
d <- d[d$subtype == "LumA" & !is.na(d$os_days) & !is.na(d$os_event) & d$os_days > 0, ]
big <- names(which(table(d$community) >= 20))
d <- d[d$community %in% big, ]
d$community <- factor(d$community)
d$nhg  <- factor(d$nhg)
d$ki67 <- factor(d$ki67_status)
cat(sprintf("LumA analizables: %d | comunidades: %s\n",
            nrow(d), paste(levels(d$community), collapse=",")))
cat("Completitud covariables: nhg no-NA =", sum(!is.na(d$nhg)),
    "| ki67 no-NA =", sum(!is.na(d$ki67)), "\n")

surv <- Surv(d$os_days, d$os_event)

# --- 2. Cox: comunidad SOLA (sin ajustar) ------------------------------------
cat("\n===== (a) Cox comunidad SOLA (sin ajuste) =====\n")
m_uni <- coxph(surv ~ community, data = d)
print(summary(m_uni)$coefficients[, c("exp(coef)","Pr(>|z|)"), drop=FALSE])

# --- 3. Cox: comunidad AJUSTADA por grado (NHG) ------------------------------
# (NHG suele ser la covariable más completa; Ki67 tiene muchos NA.)
cat("\n===== (b) Cox comunidad + NHG (grado) =====\n")
d_g <- d[!is.na(d$nhg), ]
m_grade <- coxph(Surv(os_days, os_event) ~ community + nhg, data = d_g)
print(summary(m_grade)$coefficients[, c("exp(coef)","Pr(>|z|)"), drop=FALSE])

# Test de razón de verosimilitud: ¿la comunidad aporta sobre solo NHG?
m_base <- coxph(Surv(os_days, os_event) ~ nhg, data = d_g)
lrt <- anova(m_base, m_grade)
p_lrt <- lrt$`Pr(>|Chi|)`[2]
cat(sprintf("\nLRT: ¿comunidad aporta sobre NHG solo?  p = %.4g\n", p_lrt))

# --- 4. Cox: ajustando por NHG + Ki67 (si hay datos suficientes) --------------
d_gk <- d[!is.na(d$nhg) & !is.na(d$ki67), ]
if (nrow(d_gk) >= 100 && nlevels(droplevels(d_gk$community)) >= 2) {
  cat(sprintf("\n===== (c) Cox comunidad + NHG + Ki67 (n=%d) =====\n", nrow(d_gk)))
  d_gk$community <- droplevels(d_gk$community)
  m_full <- coxph(Surv(os_days, os_event) ~ community + nhg + ki67, data = d_gk)
  print(summary(m_full)$coefficients[, c("exp(coef)","Pr(>|z|)"), drop=FALSE])
} else cat(sprintf("\n(c) omitido: pocos LumA con NHG y Ki67 (n=%d)\n", nrow(d_gk)))

# --- 5. Veredicto ------------------------------------------------------------
cat("\n===== VEREDICTO =====\n")
if (p_lrt < 0.05) {
  cat(sprintf("La comunidad topológica SIGUE siendo pronóstica tras ajustar por grado (LRT p=%.4g):\n señal independiente de los factores conocidos -> hallazgo genuino a validar.\n", p_lrt))
} else {
  cat(sprintf("La comunidad NO aporta sobre el grado (LRT p=%.4g): la topología re-deriva\n proliferación/grado -> validación de biología conocida, no novedad.\n", p_lrt))
}

save_session_info(here::here("results","sessionInfo_07_cox_luma.txt"))
