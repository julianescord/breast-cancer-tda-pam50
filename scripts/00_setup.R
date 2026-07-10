# scripts/00_setup.R
# Entorno reproducible del proyecto. Ejecutar una vez al inicio de la sesión.
# -----------------------------------------------------------------------------

# --- 1. renv: entorno de paquetes aislado y versionado -----------------------
# La PRIMERA vez en el proyecto, ejecuta manualmente en la consola:
#   install.packages("renv")
#   renv::init()          # crea renv.lock y la librería local del proyecto
# Después de instalar/actualizar paquetes, fija el estado con:
#   renv::snapshot()
# Para reproducir el entorno en otra máquina:
#   renv::restore()

# --- 2. Paquetes requeridos --------------------------------------------------
# Mirror CRAN fijo y compilación SIEMPRE desde fuente (sin binarios precompilados).
options(
  repos   = c(CRAN = "https://cloud.r-project.org"),
  pkgType = "source",
  install.packages.compile.from.source = "always"
)

# Librería personal escribible (evita necesitar sudo/root para instalar).
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && !dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

# El paquete 'fs' (>= 2.1.0) requiere libuv. Compilamos el libuv incluido para
# no depender de libuv1-dev del sistema. (Alternativa: sudo apt install libuv1-dev)
Sys.setenv(USE_BUNDLED_LIBUV = "1")

required_pkgs <- c(
  "here",          # rutas relativas a la raíz del proyecto
  "data.table",    # I/O rápido de archivos grandes
  "R.utils",       # permite a fread leer archivos .gz directamente
  "matrixStats",   # varianzas por fila eficientes
  "readr",
  "tidymodels",    # recipes, rsample, tune, workflows, yardstick, parsnip
  "kernlab",       # motor del SVM (svm_rbf / svm_poly / svm_linear)
  "FactoMineR",    # PCA (opcional; prcomp es más ligero a gran escala)
  "factoextra",
  "pheatmap"
)

installed <- rownames(installed.packages())
missing   <- setdiff(required_pkgs, installed)
if (length(missing) > 0) {
  message("Instalando (desde fuente) paquetes faltantes: ",
          paste(missing, collapse = ", "))
  install.packages(missing, type = "source")
}

invisible(lapply(required_pkgs, require, character.only = TRUE))

# --- 3. Semilla global y opciones --------------------------------------------
set.seed(2024)
options(stringsAsFactors = FALSE)

# --- 4. Funciones del proyecto -----------------------------------------------
source(here::here("R", "functions.R"))

message("Setup completo. Proyecto en: ", here::here())
