# P5 Breast Cancer — Clasificación por SVM sobre expresión génica

Estudio de machine learning para clasificar muestras de cáncer de mama a partir de
datos de expresión génica del dataset **GSE96058** (cohorte **SCAN-B**), siguiendo
una metodología estándar y reproducible.

## Estado

- [x] Estructura de proyecto reproducible
- [x] EDA + PCA (exploratorio)
- [ ] **Obtener la variable respuesta / etiquetas clínicas** (bloqueante para el SVM)
- [ ] Modelo SVM con validación cruzada (tidymodels)
- [ ] Evaluación final y reporte (DOME / TRIPOD+AI)

## Estructura

```
.
├── data/
│   ├── raw/         # datos originales (solo lectura, ignorados por git)
│   └── processed/   # derivados reproducibles
├── R/
│   └── functions.R  # funciones reutilizables
├── scripts/
│   ├── 00_setup.R       # librerías y entorno (renv)
│   ├── 01_eda_pca.R     # exploración + PCA
│   └── 02_model_svm.R   # SVM con tidymodels (sin fuga de datos)
├── results/
│   ├── figures/
│   └── tables/
├── CLAUDE.md        # metodología estándar del proyecto (leer primero)
└── README.md
```

## Cómo ejecutar

Desde **RStudio**, con el proyecto abierto en la raíz, ejecuta en orden:

```r
source("scripts/00_setup.R")     # instala/carga paquetes (renv)
source("scripts/01_eda_pca.R")   # EDA y PCA
source("scripts/02_model_svm.R") # modelo (requiere etiquetas)
```

Las rutas se resuelven con el paquete `here`, así que funciona independientemente
del directorio de trabajo mientras se ejecute dentro del proyecto.

## Datos

- **Expresión**: `data/raw/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz`
  (filas = genes, columnas = muestras; valores log2, provistos por SCAN-B).
- **Etiquetas clínicas**: aún no incluidas. Deben descargarse de GEO
  ([GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)) /
  suplementos de SCAN-B (subtipo PAM50, ER/PR/HER2, grado, supervivencia).

## Metodología

Ver [CLAUDE.md](CLAUDE.md). Puntos clave: preprocesado dentro de la validación
cruzada (sin fuga de datos), tidymodels, métricas apropiadas para clases
desbalanceadas, reproducibilidad con `renv` + semilla fija, reporte según **DOME**.
