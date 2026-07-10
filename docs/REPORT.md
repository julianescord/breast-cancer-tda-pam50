# Estratificación pronóstica del cáncer de mama Luminal A mediante análisis topológico de datos

**Proyecto P5 · GSE96058 (SCAN-B), validación en METABRIC**
Informe estructurado según el estándar **DOME** (Data, Optimization, Model, Evaluation).

---

## Resumen

Partiendo de un clasificador SVM del subtipo PAM50 (línea base reproducible pero no
novedosa), exploramos el **análisis topológico de datos (TDA)** para aportar valor.
Las *firmas topológicas por muestra* fracasaron por una razón conceptual, pero el
algoritmo **Mapper** sobre la nube de muestras (i) recuperó la biología PAM50 y
(ii) reveló un **subgrupo de Luminal A ("C2") de peor pronóstico**. Este subgrupo:

- es pronóstico **independiente de grado, proliferación e infiltración inmune** (SCAN-B);
- se caracteriza por un perfil **inmuno-frío + estromal/ECM** (COL11A1↑, marcadores T/B↓);
- **replica en una cohorte independiente** (METABRIC, otra plataforma): HR 1.23/SD, p=8×10⁻⁴;
- en METABRIC aporta sobre proliferación (≈GGI), grado y tamaño (HR 1.18, p=0.011);
- es **robusto** a los parámetros del Mapper (6/6 configuraciones) y al tamaño de la firma.

---

## D — Datos

| | SCAN-B (descubrimiento) | METABRIC (validación) |
|---|---|---|
| Fuente | GSE96058 (GEO) | cBioPortal datahub |
| Muestras | 3.409 (3.273 + 136 réplicas) | 1.980 |
| Plataforma | RNA-seq (log2) | Microarray (Illumina) |
| Genes | 30.865 | 20.384 |
| Etiquetas | PAM50 (GEO), supervivencia | PAM50, supervivencia, grado, tamaño |
| LumA analizables | 1.704 | 700 |

- **Réplicas**: `F####repl` = réplica de `F####`; agrupadas por `patient_id` para evitar fuga.
- **Etiquetas**: obtenidas de GEO con GEOquery (`03_fetch_labels.R`), mapeo 3.409/3.409.

## O — Optimización (evitar fuga de datos)

- Preprocesado que aprende de los datos (escalado, PCA) **dentro** de la validación
  cruzada (`recipes`); solo el filtro de varianza no supervisado va fuera.
- Partición **agrupada por paciente** (`group_initial_split`); test tocado una sola vez.
- Subgrupos topológicos definidos por **comunidades del grafo (sin usar supervivencia)**;
  la supervivencia se evalúa después (anti-circularidad).
- Firma de validación derivada **solo en SCAN-B** y aplicada a METABRIC sin reajuste.

## M — Modelo(s)

1. **Línea base**: SVM (RBF) para PAM50 — `02_model_svm.R`.
2. **Descubrimiento**: Mapper (Ball Mapper, `mappeR`) sobre 10 PCs → comunidades
   (Louvain) → C2 — `05_mapper.R`, `06_mapper_survival.R`.
3. **Caracterización**: expresión diferencial (limma) — `09_c2_differential_expression.R`.
4. **Confusores**: Cox ajustado por grado/proliferación/inmune — `07`,`08`,`10`.
5. **Validación**: firma portable → METABRIC — `11`,`12`.

## E — Evaluación (resultados)

### 1. Línea base SVM (validación del pipeline)
ROC-AUC 0.941, balanced accuracy 0.856 (test). Sin sobreajuste (train≈test).
Ver [svm_confusion_matrix.png](../results/figures/svm_confusion_matrix.png).
Confusión esperable LumA↔LumB; **Normal** la clase más difícil.

### 2. TDA por muestra: resultado negativo (informativo)
Firmas topológicas por muestra (imágenes de persistencia) → SVM: accuracy 0.22 (azar 0.20).
**Razón**: PAM50 se define por magnitudes de expresión de genes concretos; la homología
persistente es invariante a eso. → Descartado con fundamento.

### 3. Mapper recupera la biología PAM50
Ball Mapper: pureza media de nodos 0.74; gradiente ER (luminal→Her2→basal), flare basal.
Ver [mapper_graph.png](../results/figures/mapper_graph.png).

### 4. Subgrupo C2 (Luminal A de peor pronóstico)
Dentro de LumA, las comunidades topológicas estratifican la supervivencia
(log-rank p=1.5×10⁻⁴). C2 vs resto:

| Ajuste (Cox, SCAN-B) | HR de C2 | p |
|---|---|---|
| Sin ajustar | 2.33 | 2×10⁻⁴ |
| + grado (NHG) | 2.22 | 5×10⁻⁴ |
| + grado + proliferación | 2.16 | 1×10⁻³ |
| + grado + proliferación + inmune | 2.02 | 8×10⁻³ |

Ver [km_within_dominant.png](../results/figures/km_within_dominant.png).

### 5. Qué es C2 (expresión diferencial)
**Inmuno-frío + estromal/ECM**: bajan marcadores T (TCRB/TCRA), B/Ig (IGH/IGHV/IGJ)
y CCL19 (p≈10⁻¹⁶⁴); suben COL11A1, IBSP, HAPLN1, CEMIP.
Ver [c2_volcano.png](../results/figures/c2_volcano.png), [tabla DE](../results/tables/c2_differential_expression.csv).

### 6. Validación externa (METABRIC) y rigor
Firma C2 (SCAN-B) aplicada a 700 LumA de METABRIC → C2-like peor supervivencia:
log-rank p=0.024; **Cox continuo HR=1.23/SD, p=8×10⁻⁴**.
Ver [metabric_km_c2.png](../results/figures/metabric_km_c2.png).

| Modelo (METABRIC, HR de C2 por SD) | HR | p |
|---|---|---|
| C2 solo | 1.23 | 8×10⁻⁴ |
| + proliferación (≈GGI) | 1.23 | 1×10⁻³ |
| + grado | 1.22 | 2×10⁻³ |
| + prolif + grado + tamaño | 1.18 | 0.011 |

Robustez: firma estable con 25–200 genes (p 3×10⁻⁵ a 8×10⁻⁴).
Sensibilidad Mapper (SCAN-B): **6/6** configuraciones (eps × clustering) con log-rank
LumA significativo y comunidad peor siempre inmuno-fría. [tabla](../results/tables/scanb_sensitivity.csv).

---

## Discusión

**Qué es.** Un subgrupo de Luminal A, identificado de forma **no supervisada** por la
topología del transcriptoma, con peor pronóstico **independiente de los factores
clásicos** (grado, proliferación, inmune, tamaño) y **replicado entre plataformas**.
Su biología (inmuno-frío + estroma/ECM) es coherente con agresividad.

**Novedad y límites.** El **efecto es modesto** (HR ~1.2/SD en validación). La
contribución es doble: (a) **metodológica** — TDA no supervisada recupera estructura
pronóstica clínicamente relevante dentro de un subtipo homogéneo; (b) **una firma**
que añade sobre proliferación/grado. Para publicación faltaría: comparación formal
frente a firmas comerciales (Oncotype DX, Prosigna/ROR-P), análisis multivariable con
más covariables clínicas, y una tercera cohorte.

**Conclusión.** El hilo pasó de un ejercicio reproducible (SVM) a un **hallazgo
validado externamente**, con cada paso sometido a control (anti-fuga, anti-circular,
confusores, sensibilidad, replicación).

---

## Reproducibilidad

- Entorno: R 4.6.1, paquetes compilados desde fuente (ver `scripts/00_setup.R`); semilla `set.seed(2024)`.
- Pipeline (`scripts/`): `00` setup · `01` PCA · `03` etiquetas · `02` SVM · `04` evaluación ·
  `05`/`06` Mapper+supervivencia · `07`/`08`/`10` Cox confusores · `09` DE ·
  `11`/`12` validación METABRIC · `13` sensibilidad.
- Salidas en `results/` (figuras, tablas, `sessionInfo_*`). Datos crudos en `data/raw/` (no versionados).
