# P5 Breast Cancer — Guía del proyecto

Estudio de **machine learning predictivo** (Support Vector Machine) sobre datos de
expresión génica de cáncer de mama (**GSE96058 / cohorte SCAN-B**). El objetivo es
un flujo **estándar y reproducible**, no un pipeline de expresión diferencial clásico.

Este archivo define cómo se trabaja en el proyecto. Cualquier análisis nuevo debe
respetar estas reglas.

---

## Entorno

- El código es **R**, ejecutado desde **RStudio** (no hay R en el PATH de la terminal;
  no asumas que puedes correr `Rscript` — los scripts se ejecutan manualmente en RStudio).
- Rutas **siempre relativas a la raíz del proyecto** usando el paquete `here`
  (`here::here("data","raw", ...)`). Nunca rutas absolutas ni `setwd()`.
- Gestión de paquetes con **`renv`** (ver `scripts/00_setup.R`). Tras instalar o
  actualizar un paquete: `renv::snapshot()` para fijarlo en `renv.lock`.
- **Semilla fija** en todo lo estocástico: `set.seed(2024)` (o la definida en el script).
- Cada análisis termina registrando `sessionInfo()` en `results/`.

## Estructura de carpetas

| Carpeta            | Contenido                                                        |
|--------------------|------------------------------------------------------------------|
| `data/raw/`        | Datos originales, **solo lectura**. Nunca se sobrescriben.        |
| `data/processed/`  | Derivados reproducibles (matrices filtradas, PCA, splits).       |
| `R/`               | Funciones reutilizables (`functions.R`), sin efectos secundarios.|
| `scripts/`         | Scripts numerados que se ejecutan en orden (`00_`, `01_`, ...).  |
| `results/figures/` | Figuras exportadas.                                              |
| `results/tables/`  | Tablas y métricas exportadas.                                    |

## Datos

- **Matriz de expresión**: `data/raw/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz`
  - Orientación: **filas = genes**, **columnas = muestras** (`F1`…). Primera columna = ID de gen.
  - Para modelar hay que **transponer**: observaciones = muestras, variables = genes.
  - Valores ya transformados (log2) por los autores de SCAN-B.
- **Etiquetas / desenlace**: se obtienen de GEO con `scripts/03_fetch_labels.R`
  (GEOquery), que genera `data/raw/metadata.csv` mapeando por `sample` (F####,
  = `title` de GEO). Variable respuesta del estudio: **subtipo PAM50** (`pam50`,
  multiclase: LumA/LumB/Her2/Basal/Normal; desbalanceado, LumA ~50%). El CSV
  incluye también er/pgr/her2/ki67, nhg y supervivencia para análisis futuros.
- **Réplicas**: hay 136 réplicas técnicas (`F####repl` = réplica de `F####`). Del
  mismo tumor **no pueden** repartirse entre train/test (fuga). Se agrupan por
  `patient_id` (= `F####` sin sufijo) en el split (`group_initial_split` / `group_vfold_cv`).

---

## Pipeline (scripts en orden)

`00_setup` (entorno) → `01_eda_pca` (PCA exploratorio) → `03_fetch_labels`
(etiquetas PAM50 desde GEO) → `02_model_svm` (entrena+tunea SVM) → `04_evaluate`
(métricas por clase train/test).

```
DATOS  GSE96058: 30.865 genes × 3.409 muestras (log2)  +  PAM50 (metadata.csv)
  │
1·PREPROCESO (fuera de CV, NO supervisado, sin etiqueta):
    filtro varianza 30.865→3.000 genes · transponer · unir PAM50 · split 80/20 group_by(paciente)
  ══════ frontera anti-fuga: todo lo de abajo se reajusta en el train de cada fold ══════
2·RECETA (en CV):  step_zv → step_normalize → step_pca(num_comp)
3·MODELO:  SVM kernel RBF · tuning CV 5-fold × 15 combos · mejor por ROC-AUC · fit final en train
4·RESULTADOS (test 1 vez):  ROC-AUC, accuracy, bal_acc, kappa · por clase: prec/recall/spec/F1 · confusión
```
Nota: el SVM **no tiene epochs** (optimización convexa de una sola solución); el
análogo a "curva de entrenamiento" es la búsqueda de hiperparámetros por CV
(`results/figures/svm_tuning.png`, `results/tables/svm_tuning_metrics.csv`). ANOVA
NO es el análisis de resultados de una clasificación (compara medias entre grupos);
solo aplicaría como selección supervisada de genes, y entonces iría DENTRO de la CV.

## Resultados (SVM PAM50 — corrida completa: 3000 genes, 5 folds, 15 combos)

Mejores hiperparámetros: `cost=21.5`, `rbf_sigma=2e-4`, `num_comp=39`.

Global (train ≈ test → **sin sobreajuste**):

| set   | accuracy | macro-F1 | kappa | ROC-AUC | recall medio |
|-------|----------|----------|-------|---------|--------------|
| train | 0.855    | 0.812    | 0.777 | 0.952   | 0.780        |
| test  | 0.850    | 0.795    | 0.769 | 0.941   | 0.764        |

Por clase (TEST):

| clase  | n   | precisión | recall/sensib. | especif. | F1    |
|--------|-----|-----------|----------------|----------|-------|
| Basal  | 86  | 0.975     | 0.919          | 0.997    | 0.946 |
| Her2   | 60  | 0.864     | 0.850          | 0.987    | 0.857 |
| LumA   | 341 | 0.835     | 0.933          | 0.815    | 0.881 |
| LumB   | 156 | 0.818     | 0.750          | 0.951    | 0.783 |
| Normal | 38  | 0.824     | **0.368**      | 0.995    | 0.509 |

Matriz de confusión TEST (filas=predicción, cols=verdad):

```
        Basal Her2 LumA LumB Normal
Basal      79    0    0    1     1
Her2        2   51    3    1     2
LumA        2    3  318   37    21
LumB        1    6   19  117     0
Normal      2    0    1    0    14
```
Lectura: Basal/Her2 muy bien; confusión esperable LumA↔LumB; **Normal** es la clase
difícil (recall 0.37, se confunde con LumA). Aviso: `yardstick::bal_accuracy`=0.856
usa (sens+espec)/2 (la alta especificidad la infla); el "recall medio"=0.764 (media
de sensibilidad) es la lectura estricta. Fuentes: `results/tables/svm_overall_train_test.csv`,
`svm_perclass_{train,test}.csv`, `svm_confusion_{train,test}.csv`.

---

## Reglas metodológicas (obligatorias)

### 1. Evitar fuga de datos (data leakage) — la regla más importante
Todo preprocesado que **aprende parámetros de los datos** (escalado, PCA, selección
de genes por varianza/estadístico, imputación, balanceo) debe ajustarse **solo con
el conjunto de entrenamiento dentro de cada fold** y aplicarse al test.
- **NO** hacer PCA/escalado/selección sobre la matriz completa antes de dividir.
- En la práctica: implementarlo como pasos de una `recipe` de **tidymodels**
  (`step_normalize`, `step_pca`, `step_zv`, `step_nzv`, …), que `tidymodels`
  reajusta automáticamente en cada remuestreo.
- **Excepción aceptada**: un filtro de genes por varianza **no supervisado** (que
  NO usa la etiqueta) puede hacerse fuera de la CV; es práctica estándar y necesario
  por coste computacional (30k genes → ~3000). Cualquier filtro que use el desenlace
  (t-test, correlación con la clase, etc.) SÍ debe ir dentro de la CV.

### 2. Partición y validación
- `initial_split` (o `group_initial_split` por paciente) → train/test, estratificado
  por la clase objetivo.
- Ajuste de hiperparámetros del SVM (coste `C`, `sigma` en RBF) con **validación
  cruzada** (`vfold_cv`, idealmente estratificada); para estimaciones honestas
  considerar **CV anidada** (`nested_cv`).
- El **test set se toca una sola vez**, al final (`last_fit`).

### 3. Métricas
- No usar solo *accuracy* (clases desbalanceadas). Reportar AUC-ROC, sensibilidad,
  especificidad, F1, matriz de confusión, y balanced accuracy (`yardstick`).

### 4. Reproducibilidad
- Semilla fija, `renv.lock` versionado, scripts numerados y deterministas,
  `sessionInfo()` guardado, todo bajo control de versiones (git).

### 5. Reporte según estándares
- Documentar el estudio siguiendo **DOME** (Data, Optimization, Model, Evaluation;
  Nat. Methods 2021) y, si el desenlace es clínico, **TRIPOD+AI**.

---

## Convenciones de código
- Estilo tidyverse; nombres en `snake_case`.
- Funciones puras en `R/functions.R`; los `scripts/` orquestan y guardan salidas.
- Preferir `data.table::fread`/`readr` para I/O; para PCA a gran escala usar
  `prcomp` sobre genes de alta varianza (FactoMineR es costoso con ~30k genes).
