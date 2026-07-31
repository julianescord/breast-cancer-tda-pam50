# scripts/03_fetch_labels.R
# Descarga los metadatos clínicos de GSE96058 (SCAN-B) desde GEO y construye
# data/raw/metadata.csv con la variable respuesta (subtipo PAM50) y campos
# clínicos auxiliares, mapeando por el identificador de muestra (F####).
#
# Requiere: source("scripts/00_setup.R") y GEOquery (Bioconductor).
# -----------------------------------------------------------------------------

source(here::here("scripts", "00_setup.R"))
suppressMessages(library(GEOquery))

# Caché de la descarga dentro de data/raw (ignorada por git).
geo_cache <- here::here("data", "raw", "geo_cache")
dir.create(geo_cache, recursive = TRUE, showWarnings = FALSE)

# --- 1. Descargar metadatos (2 plataformas: HiSeq 2000 y NextSeq 500) --------
gse <- getGEO("GSE96058", GSEMatrix = TRUE, getGPL = FALSE, destdir = geo_cache)

# --- 2. Extraer un campo de 'characteristics_ch1*' por su etiqueta -----------
# Los campos vienen como "etiqueta: valor". Buscamos por etiqueta (robusto al
# orden/índice de columnas, que puede variar entre plataformas).
get_char <- function(pd, label) {
  cc <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
  pat <- paste0("^", label, ":")
  for (col in cc) {
    vals <- as.character(pd[[col]])
    if (any(grepl(pat, vals))) {
      out <- sub(paste0("^", label, ":\\s*"), "", vals)
      out[!grepl(pat, vals)] <- NA  # por si algún valor no corresponde
      return(trimws(out))
    }
  }
  rep(NA_character_, nrow(pd))
}

# --- 3. Construir tabla de etiquetas por plataforma --------------------------
build_meta <- function(eset) {
  pd <- Biobase::pData(eset)
  data.frame(
    sample       = trimws(as.character(pd$title)),
    pam50        = get_char(pd, "pam50 subtype"),
    er_status    = get_char(pd, "er status"),
    pgr_status   = get_char(pd, "pgr status"),
    her2_status  = get_char(pd, "her2 status"),
    ki67_status  = get_char(pd, "ki67 status"),
    nhg          = get_char(pd, "nhg"),
    os_days      = get_char(pd, "overall survival days"),
    os_event     = get_char(pd, "overall survival event"),
    instrument   = get_char(pd, "instrument model"),
    stringsAsFactors = FALSE
  )
}

meta <- do.call(rbind, lapply(gse, build_meta))

# --- 4. patient_id para agrupar réplicas (F3028 y F3028repl = misma muestra) --
meta$patient_id <- sub("repl$", "", meta$sample)

# --- 5. Verificación contra las columnas de la matriz de expresión -----------
expr_path <- here::here(
  "data", "raw",
  "GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
)
expr_cols <- names(data.table::fread(expr_path, nrows = 0))[-1]  # sin col de genes

matched <- intersect(expr_cols, meta$sample)
message(sprintf("Muestras en expresión: %d | en metadatos: %d | coinciden: %d",
                length(expr_cols), nrow(meta), length(matched)))
missing_expr <- setdiff(expr_cols, meta$sample)
if (length(missing_expr) > 0) {
  warning("Columnas de expresión SIN metadato: ",
          paste(head(missing_expr, 10), collapse = ", "),
          if (length(missing_expr) > 10) " ..." else "")
}

# --- 6. Distribución de la clase objetivo ------------------------------------
cat("\nDistribución de PAM50:\n")
print(table(meta$pam50, useNA = "ifany"))

# --- 7. Guardar --------------------------------------------------------------
readr::write_csv(meta, here::here("data", "raw", "metadata.csv"))
message("\nmetadata.csv escrito en data/raw/. Variable respuesta: `pam50`.")
save_session_info(here::here("results", "sessionInfo_03_fetch_labels.txt"))
