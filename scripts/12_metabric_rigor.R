# scripts/12_metabric_rigor.R
# RIGOR en METABRIC: ¿la firma C2 sigue siendo pronóstica en LumA tras ajustar por
# GRADO, PROLIFERACIÓN (~GGI) y TAMAÑO?  ¿Aporta sobre la proliferación (la firma
# luminal pronóstica clásica)?  ¿Es robusta al tamaño de la firma?
#
# Requiere: 09 (DE), descargas METABRIC (expresión + clínico paciente y muestra).
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(survival))
ext <- here::here("data","raw","external")

# --- 1. Clínico: subtipo + supervivencia + grado + tamaño --------------------
clp <- data.table::fread(file.path(ext,"metabric_clinical_patient.txt"), skip="PATIENT_ID")[
         , c("PATIENT_ID","CLAUDIN_SUBTYPE","OS_MONTHS","OS_STATUS")]
cls <- data.table::fread(file.path(ext,"metabric_clinical_sample.txt"), skip="PATIENT_ID")[
         , c("PATIENT_ID","GRADE","TUMOR_SIZE")]
cl <- merge(clp, cls, by="PATIENT_ID")
cl$event    <- as.integer(grepl("DECEASED|^1", cl$OS_STATUS))
cl$OS_MONTHS <- suppressWarnings(as.numeric(cl$OS_MONTHS))
cl$GRADE     <- factor(cl$GRADE)
cl$TUMOR_SIZE <- suppressWarnings(as.numeric(cl$TUMOR_SIZE))

# --- 2. Expresión METABRIC ---------------------------------------------------
ex <- data.table::fread(file.path(ext,"metabric_expression.txt"))
sym <- ex[[1]]; Xm <- as.matrix(ex[, -(1:2)])
Xm <- rowsum(Xm, sym, na.rm=TRUE) / as.numeric(table(sym)[rownames(rowsum(Xm, sym))])
Zm <- t(scale(t(Xm)))

# --- 3. Scores: firma C2 (de SCAN-B) y proliferación -------------------------
de <- readr::read_csv(here::here("results","tables","c2_differential_expression.csv"),
                      show_col_types=FALSE)
de <- de[de$adj.P.Val < 0.05, ]
mb <- ifelse(de$gene %in% rownames(Xm), de$gene,
             ifelse(gsub(".", "-", de$gene, fixed=TRUE) %in% rownames(Xm),
                    gsub(".", "-", de$gene, fixed=TRUE), NA))
de$mb <- mb; de <- de[!is.na(de$mb), ]
c2_sig <- function(n) {
  up   <- head(de$mb[order(-de$logFC)], n); down <- head(de$mb[order(de$logFC)], n)
  as.numeric(scale(colMeans(Zm[up,,drop=FALSE],na.rm=TRUE) - colMeans(Zm[down,,drop=FALSE],na.rm=TRUE)))
}
prolif_g <- intersect(c("ANLN","AURKA","BIRC5","CCNB1","CCNB2","CDC20","CDK1","CDKN3",
  "CENPF","CEP55","EXO1","FOXM1","KIF2C","MELK","MKI67","NDC80","NUF2","PTTG1","RRM2",
  "TOP2A","TYMS","UBE2C"), rownames(Xm))
prolif <- as.numeric(scale(colMeans(Zm[prolif_g,,drop=FALSE], na.rm=TRUE)))

sc <- data.frame(PATIENT_ID=colnames(Xm), c2=c2_sig(100), prolif=prolif)
d <- merge(cl, sc, by="PATIENT_ID")
d <- d[d$CLAUDIN_SUBTYPE=="LumA" & !is.na(d$OS_MONTHS) & d$OS_MONTHS>0 & !is.na(d$event), ]
cat(sprintf("LumA METABRIC: %d | con grado no-NA: %d\n", nrow(d), sum(!is.na(d$GRADE))))

# --- 4. Modelos de Cox (HR por SD del score C2) ------------------------------
show_c2 <- function(m) { s <- summary(m)$coefficients; s["c2", c("exp(coef)","Pr(>|z|)")] }
cat("\n=== C2 (HR por SD) en modelos crecientes ===\n")
cat(sprintf("(a) C2 solo:              HR=%.2f  p=%.4g\n", show_c2(coxph(Surv(OS_MONTHS,event)~c2, d))[1], show_c2(coxph(Surv(OS_MONTHS,event)~c2, d))[2]))
m_p <- coxph(Surv(OS_MONTHS,event)~c2+prolif, d)
cat(sprintf("(b) C2 + proliferación:   HR=%.2f  p=%.4g   [prolif p=%.3g]\n", show_c2(m_p)[1], show_c2(m_p)[2], summary(m_p)$coefficients["prolif","Pr(>|z|)"]))
dg <- d[!is.na(d$GRADE), ]
m_g <- coxph(Surv(OS_MONTHS,event)~c2+GRADE, dg)
cat(sprintf("(c) C2 + grado:           HR=%.2f  p=%.4g\n", show_c2(m_g)[1], show_c2(m_g)[2]))
dgs <- d[!is.na(d$GRADE) & !is.na(d$TUMOR_SIZE), ]
m_full <- coxph(Surv(OS_MONTHS,event)~c2+prolif+GRADE+TUMOR_SIZE, dgs)
ci_full <- summary(m_full)$conf.int["c2", c("lower .95","upper .95")]
zph_full <- as.numeric(survival::cox.zph(m_full)$table["GLOBAL","p"])
cat(sprintf("(d) C2 + prolif+grado+tam:HR=%.2f (IC95%% %.2f-%.2f)  p=%.4g   (n=%d)  [PH global p=%.3g]\n",
            show_c2(m_full)[1], ci_full[1], ci_full[2], show_c2(m_full)[2], nrow(dgs), zph_full))

# --- 5. Sensibilidad al tamaño de la firma -----------------------------------
cat("\n=== Sensibilidad: HR/p de C2 según nº de genes por lado ===\n")
for (n in c(25,50,100,200)) {
  d$c2n <- c2_sig(n)[match(d$PATIENT_ID, colnames(Xm))]
  m <- coxph(Surv(OS_MONTHS,event)~c2n, d)
  s <- summary(m)$coefficients
  cat(sprintf("  %3d genes/lado: HR=%.2f  p=%.4g\n", n, s["c2n","exp(coef)"], s["c2n","Pr(>|z|)"]))
}
save_session_info(here::here("results","sessionInfo_12_metabric_rigor.txt"))
