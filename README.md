# P5 Breast Cancer — Clasificación PAM50 y descubrimiento de subgrupo pronóstico por TDA

Clasificación PAM50 con SVM y descubrimiento de un subgrupo pronóstico de Luminal A
mediante Mapper (TDA) a partir de datos de expresión génica del dataset **GSE96058**
(cohorte **SCAN-B**), validado externamente en **METABRIC**. Pipeline reproducible
en R/tidymodels.

## Estado

- [x] Estructura de proyecto reproducible
- [x] EDA + PCA (exploratorio)
- [x] Etiquetas clínicas (subtipo PAM50, ER/PR/HER2, grado, supervivencia)
- [x] Modelo SVM con validación cruzada (tidymodels): **85.0% accuracy, ROC AUC 0.94
      (test set)**, sin fuga de datos
- [x] Evaluación final (matriz de confusión, métricas por clase)
- [x] Descubrimiento de subgrupos vía Mapper (TDA) y análisis de supervivencia (Cox)
- [x] Expresión diferencial y análisis de inmunidad/proliferación por subgrupo
- [x] Validación externa en METABRIC
- [x] Análisis de sensibilidad
- [ ] Reporte final (DOME / TRIPOD+AI)

## Estructura

```
.
├── data/
│   ├── raw/         # datos originales (solo lectura, ignorados por git)
│   └── processed/   # derivados reproducibles
├── R/
│   └── functions.R  # funciones reutilizables
├── scripts/
│   ├── 00_setup.R                      # librerías y entorno (renv)
│   ├── 01_eda_pca.R                    # exploración + PCA
│   ├── 02_model_svm.R                  # SVM con tidymodels (sin fuga de datos)
│   ├── 03_fetch_labels.R               # etiquetas clínicas (PAM50, ER/PR/HER2, supervivencia)
│   ├── 04_evaluate.R                   # evaluación del modelo (test set)
│   ├── 05_mapper.R                     # descubrimiento de subgrupos (TDA/Mapper)
│   ├── 06_mapper_survival.R            # supervivencia por comunidad Mapper
│   ├── 07_cox_luma.R                   # modelo de Cox en el subgrupo Luminal A
│   ├── 08_proliferation_cox.R          # Cox ajustado por proliferación
│   ├── 09_c2_differential_expression.R # expresión diferencial del subgrupo pronóstico
│   ├── 10_immune_cox.R                 # Cox ajustado por infiltración inmune
│   ├── 11_external_validation.R        # validación en METABRIC
│   ├── 12_metabric_rigor.R             # chequeos de rigor en METABRIC
│   ├── 13_sensitivity_scanb.R          # análisis de sensibilidad
│   ├── 14_build_report_html.R          # reporte HTML
│   └── 15_figures_english.R            # figuras en inglés para publicación
├── results/
│   ├── figures/
│   └── tables/
├── CLAUDE.md        # metodología estándar del proyecto (leer primero)
└── README.md
```

## Cómo ejecutar

Desde **RStudio**, con el proyecto abierto en la raíz, ejecuta los scripts en
orden numérico (`00_setup.R` a `15_figures_english.R`). Las rutas se resuelven
con el paquete `here`, así que funciona independientemente del directorio de
trabajo mientras se ejecute dentro del proyecto.

## Datos

- **Expresión**: `data/raw/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz`
  (filas = genes, columnas = muestras; valores log2, provistos por SCAN-B).
- **Etiquetas clínicas**: subtipo PAM50, ER/PR/HER2, grado y supervivencia,
  obtenidas de GEO ([GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058))
  y suplementos de SCAN-B.
- **Validación externa**: cohorte METABRIC.

## Metodología

Ver [CLAUDE.md](CLAUDE.md). Puntos clave: preprocesado dentro de la validación
cruzada (sin fuga de datos), tidymodels, métricas apropiadas para clases
desbalanceadas, reproducibilidad con `renv` + semilla fija, reporte según **DOME**.
