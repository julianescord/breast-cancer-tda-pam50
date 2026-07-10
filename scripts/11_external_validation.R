# scripts/11_external_validation.R
# VALIDACIÓN EXTERNA en METABRIC: ¿los LumA "C2-like" (firma portable derivada de
# la DE de SCAN-B) también tienen peor supervivencia en una cohorte independiente?
#
# Firma C2 = media(z de genes SUBEN en C2) - media(z de genes BAJAN en C2).
# Se deriva en SCAN-B, se aplica a METABRIC (nada de esta cohorte se usó antes).
#
# Requiere: 09_c2_differential_expression.R (DE) + descarga METABRIC. Paquete: survival.
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(survival))
ext <- here::here("data","raw","external")

# --- 1. Clínico METABRIC: subtipo PAM50 + supervivencia ----------------------
cl <- data.table::fread(file.path(ext,"metabric_clinical_patient.txt"), skip="PATIENT_ID")
cl <- cl[, c("PATIENT_ID","CLAUDIN_SUBTYPE","OS_MONTHS","OS_STATUS")]
cl$event <- as.integer(grepl("DECEASED|^1", cl$OS_STATUS))
cl$OS_MONTHS <- suppressWarnings(as.numeric(cl$OS_MONTHS))
cat("Subtipos METABRIC:\n"); print(table(cl$CLAUDIN_SUBTYPE))

# --- 2. Expresión METABRIC (genes x muestras) --------------------------------
ex <- data.table::fread(file.path(ext,"metabric_expression.txt"))
sym <- ex[[1]]                                    # Hugo_Symbol
Xm <- as.matrix(ex[, -(1:2)])                     # quita Hugo + Entrez
# colapsar símbolos duplicados por media
Xm <- rowsum(Xm, sym, na.rm=TRUE) / as.numeric(table(sym)[rownames(rowsum(Xm, sym))])
message(sprintf("Expresión METABRIC: %d genes x %d muestras", nrow(Xm), ncol(Xm)))

# --- 3. Firma portable de C2 (desde la DE de SCAN-B) -------------------------
de <- readr::read_csv(here::here("results","tables","c2_differential_expression.csv"),
                      show_col_types = FALSE)
de <- de[de$adj.P.Val < 0.05, ]
# nombres de SCAN-B pasaron por make.names (guiones -> puntos); probamos ambos
to_metabric <- function(g) { g2 <- gsub(".", "-", g, fixed=TRUE)
  ifelse(g %in% rownames(Xm), g, ifelse(g2 %in% rownames(Xm), g2, NA)) }
de$mb <- to_metabric(de$gene); de <- de[!is.na(de$mb), ]
up   <- head(de$mb[order(-de$logFC)], 100)        # suben en C2
down <- head(de$mb[order( de$logFC)], 100)        # bajan en C2
cat(sprintf("Genes de la firma presentes en METABRIC: %d up, %d down\n", length(up), length(down)))

Zm <- t(scale(t(Xm)))                              # z-score por gen
c2_score <- colMeans(Zm[up, , drop=FALSE], na.rm=TRUE) - colMeans(Zm[down, , drop=FALSE], na.rm=TRUE)
sig <- data.frame(PATIENT_ID = colnames(Xm), c2_score = as.numeric(scale(c2_score)))

# --- 4. Unir y quedarnos con LumA --------------------------------------------
d <- merge(cl, sig, by = "PATIENT_ID")
d <- d[d$CLAUDIN_SUBTYPE == "LumA" & !is.na(d$OS_MONTHS) & !is.na(d$event) & d$OS_MONTHS > 0, ]
d$c2like <- factor(ifelse(d$c2_score > median(d$c2_score), "C2-like", "resto"),
                   levels = c("resto","C2-like"))
cat(sprintf("\nLumA en METABRIC analizables: %d (C2-like=%d, resto=%d)\n",
            nrow(d), sum(d$c2like=="C2-like"), sum(d$c2like=="resto")))

# --- 5. Supervivencia: ¿C2-like peor en METABRIC? ----------------------------
lr <- survdiff(Surv(OS_MONTHS, event) ~ c2like, data = d)
p_lr <- 1 - pchisq(lr$chisq, 1)
cox <- coxph(Surv(OS_MONTHS, event) ~ c2_score, data = d)          # continuo, por SD
hr  <- summary(cox)$coefficients[, "exp(coef)"]
p_cox <- summary(cox)$coefficients[, "Pr(>|z|)"]
cat(sprintf("\nLog-rank C2-like vs resto (LumA METABRIC): p = %.4g\n", p_lr))
cat(sprintf("Cox por score continuo: HR = %.2f por SD, p = %.4g\n", hr, p_cox))

png(here::here("results","figures","metabric_km_c2.png"), width=1100, height=800, res=150)
plot(survfit(Surv(OS_MONTHS, event) ~ c2like, data=d), col=c("#4575b4","#d73027"), lwd=2,
     xlab="meses", ylab="Supervivencia",
     main=sprintf("METABRIC · LumA C2-like vs resto (log-rank p=%.3g)", p_lr))
legend("bottomleft", legend=levels(d$c2like), col=c("#4575b4","#d73027"), lwd=2, bty="n")
dev.off()

cat("\n===== VEREDICTO VALIDACIÓN EXTERNA =====\n")
if (p_lr < 0.05) cat(sprintf("REPLICA: los LumA C2-like también tienen peor supervivencia en METABRIC (p=%.4g).\n", p_lr)) else cat(sprintf("NO replica en METABRIC (p=%.4g): el hallazgo puede ser específico de SCAN-B.\n", p_lr))
save_session_info(here::here("results","sessionInfo_11_external_validation.txt"))
