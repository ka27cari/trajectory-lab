# Trajectory Lab 🧬

**An interactive Shiny sandbox for synthetic methylation-like time courses, stepwise PCA, controlled mean/variance changes, feature correlations, and distance-based embedding.**

Generate bounded (0–1) logistic trajectories, change the observation period, inspect evolving PCA scores and loadings, compare four 2D maps at once, click a correlation heatmap, upload your own beta-scale matrix, and export the results. This is a standalone rewrite inspired by the ideas in two older R/Shiny research prototypes, not a redistribution of their raw research data.

> **Scientific scope:** The simulated curves are illustrative, not biological observations or validated methylation estimates. The app does not perform EWAS, batch correction, cell-type deconvolution, causal inference, or clinical prediction.

## Features

| Area | Capabilities |
| --- | --- |
| Simulator | Mixed, early, late, and synchronized transition patterns; amplitude, timing, noise, feature count, and random seed controls |
| Editor | Add custom sigmoid trajectories and undo the last addition |
| Import | Validate and analyze an uploaded beta-scale CSV in session memory |
| Explorer | Highlight any feature and inspect its trajectory and generator metadata |
| PCA | Standardized PCA, score trajectories, explained variance, PC1 loadings, **feature-by-feature PCA refit with orthogonal alignment**, fixed-basis partial contributions, and a cumulative loading walk |
| Mean & variance | Compare earlier/later sample halves by mean change and variance ratio; independently control a feature’s mean offset and spread, optionally applying it to PCA, distance maps and correlation |
| Distance maps | **Four views together**: PCA (original values), PCoA, nonmetric MDS, and Sammon mapping; five pairwise distances; map either time points or features; quality diagnostics, matrix heatmap and 2D distance plots |
| Correlations | Pearson or Spearman; optional Gaussian time weights; clickable heatmap; pairwise scatter and time-course views |
| Export | Original and active-window data, metadata, PCA and stepwise coordinates, correlation matrix, distance matrix, four-method coordinates and quality table, settings, and PNG plots |
| Quality | Pure analysis functions, reproducibility tests, GitHub Actions workflow, explicit assumptions |

**Version history:** This v1.2.0 source continues the previously interrupted v1.1.0 work and retains the v1.0.1 correlation-matrix shape fix (`MEAN |r| = N/A` / `dimnames` errors). Restart Shiny after replacing files. Calculation and chart code have been statically inspected, **not run locally**, because R is unavailable in this environment. GitHub Actions is configured to run the tests after publishing.

## Quick start (R)

Requires a recent R installation, an internet connection when installing packages, and the packages **shiny**, **bslib**, **ggplot2**, and **MASS** (a recommended R package used for NMDS/Sammon).

```r
# From the extracted project folder in R or RStudio:
source("install.R")
shiny::runApp(".")
```

To verify the independent calculation functions:

```sh
Rscript tests/run_tests.R
```

Pure computational tests use base R and are automatically run by `.github/workflows/r-tests.yml` after pushing to GitHub. CI also installs the UI and ordination dependencies and verifies that `source("app.R")` constructs a Shiny app. This is **not** a browser end-to-end test; manually check chart rendering and responsiveness using the Korean smoke-test checklist.

## Input CSV format

The first file row must contain headers. A numeric column exactly named `time` is required, along with 2–40 columns containing beta-scale feature values between `0` and `1`.

```csv
time,CpG_01,CpG_02,CpG_03
0,0.42,0.73,0.29
100,0.44,0.71,0.31
200,0.46,0.69,0.32
300,0.48,0.65,0.34
400,0.50,0.62,0.36
```

Try the included [`examples/demo-beta.csv`](examples/demo-beta.csv) to exercise the import tab. This demonstration file contains only synthetic values.

Import constraints: **5–1000 rows**, **2–40 features**, **3 MB maximum**, strictly increasing time, finite numeric values, no missing cells, and no duplicate column names. Invalid data are rejected rather than silently imputed. Uploaded datasets are read into the current session. The application code does not persist them to a database or write copies; inspect your hosting provider's handling and retention policies before uploading sensitive data.

## Feature-by-feature PCA and shift/spread experiments (v1.1.0)

**PCA → Feature-by-feature PCA.** Move the **Features included** slider (or click **Previous / Add next feature**) from 2 features to the entire feature set in refit mode, or from 1 feature in fixed-basis mode. Pick the inclusion order: original CSV/creation order, late-to-early transition midpoint, or descending absolute PC1 loading. **Adding a third principal component to an already fitted PCA does not move its PC1/PC2 scores.** This stepper adds *features/variables* to a fit, or their contributions to fixed axes; it does not pretend that including PC3 changes the first two axes. The two modes answer different questions:

- **Refit and align PCA each step:** recomputes `prcomp()` from scratch with the selected features. The new 2D score map is orthogonally aligned to the full-feature 2D reference (Procrustes alignment) for legibility. **Different fits do not have literally identical PC axes**, and the scores are not directly equivalent to a fixed projection. Grey shows the full-feature reference.
- **Accumulate contributions in fixed axes:** computes partial scores by adding standardized feature × its loading in a fixed full-data PCA basis. Only the contribution grows, not the fitted principal directions. This mode starts at one feature; a PCA refit begins with two features.
- **Cumulative loading walk:** starts at `(0, 0)` and adds the PC1/PC2 loading vector of each selected feature in order. It is a trajectory **of loading-vector sums, not time-point score locations**. This preserves the visually exploratory intent of the older prototype without claiming both quantities are the same.

**Mean & variance tab.** The selected observation window is split by **sample order into two consecutive halves** (not independent subjects). The summary compares `late mean − early mean` and `late variance / early variance` for each feature; a zero early variance yields `NA` rather than a misleading ratio. A variance ratio below 1 means narrower variation in the later half; means and variances can change independently.

In the transformation lab, pick a feature and try **Requested mean offset** with **Spread multiplier**. For baseline feature values `x`, the transform is `clamp(mean(x) + offset + multiplier × (x − mean(x)), 0, 1)`. A mean-only shift does not change variance *unless clipping occurs*; a multiplier of 0.5 without clipping yields 0.25× variance. The app reports the **realized** shift, variance ratio, and number clipped. Changes stay in preview until **Apply to PCA, distance maps, correlations and exports** is enabled. The original full dataset and metadata exports are preserved; active-window, PCA, and correlation exports reflect the enabled experiment.

**Chart legibility:** continuous legends use a compact horizontal color bar with sparse, shortened ticks. The correlation heatmap selectively displays x/y labels when many features are present, preserving the entire matrix and click targets rather than drawing overlapping text.

## Distance maps: beyond PCA (v1.2.0)

In the **Distance maps** tab, set a pairwise distance, select which objects to map, and press **Recalculate four maps**. The app displays **PCA / PCoA / NMDS / Sammon side-by-side**, then a common diagnostic table, detailed map, a pairwise distance-preservation plot, and a thinned, readable distance heatmap.

| Control | What changes |
| --- | --- |
| **Distance** | Euclidean (L2), Manhattan (L1), Canberra, Bray–Curtis, or `1 - Pearson r` between object profiles. **PCA does not use this choice**; PCoA, NMDS, and Sammon do. |
| **Objects** | Rows = time points across features (default), or rows = features across sampled times. The latter maps relationships between entire feature profiles, not trajectories of time points. |
| **Standardize** | Center and z-score each nonconstant feature across time before computing distances/embeddings. Turn **off** for Bray–Curtis; signed standardized values invalidate it. The feature-orientation option uses the same per-feature preprocessing first, then transposes. |
| **Point cap** | For responsive server-side MDS, when mapping time points the app **evenly samples at most 40–200 time points** (default 120) from the active analysis window. The remaining tabs retain their full active data. Feature mapping retains all eligible feature profiles (at most 40). |
| **Recalculate** | An explicit trigger for potentially expensive MDS fits. Any data, experiment, window, distance or preprocess changes are incorporated on the next click; the charts and exports otherwise retain the last computed maps. |

Methods:

- **PCA**: center the prepared objects and extract principal component scores; compares an original-value matrix to the distance maps. Changes in *distance metric alone* do not alter the PCA coordinates.
- **PCoA / classical MDS**: `stats::cmdscale()` finds 2D coordinates from the selected distance matrix. For non-Euclidean distances, negative eigenvalues may occur; the method detail reports the negative eigenvalue mass and positive-inertia fraction. These are *not* PCA explained variances.
- **NMDS**: `MASS::isoMDS()` approximately preserves the **rank ordering** of the dissimilarities through a monotone fit (not necessarily their magnitudes). Reports the method's own Kruskal stress, in percent.
- **Sammon**: `MASS::sammon()` puts extra weight on small dissimilarities in the optimization and reports its own Sammon stress. It is a nonlinear distance-based map.

The overview table displays **distance-rank correlation** (Spearman's rho between the selected dissimilarities and 2D Euclidean distances, higher is better) and **relative distance error** (residual after fitting a single nonnegative global scale, lower is better). These are *descriptive and not equivalent to each algorithm's objective*. In particular, a strong NMDS rank score can coexist with a relatively large metric distance error. The diagnostic scatter uses the globally rescaled distances. 2D maps may rotate or reflect, so compare their local structure and not raw `(x, y)` coordinates.

Input rows with identical profiles produce zero off-diagonal dissimilarities, which `isoMDS()` and `sammon()` do not accept. Those maps display a specific **Unavailable** message without preventing PCA/PCoA display. A constant profile makes correlation distance undefined; an all-zero pair makes Bray–Curtis undefined. The app reports these instead of inventing output or jittering observations. Each method uses a deterministic classical MDS initialization, though optimization results can still vary across R/library versions.

**Important:** This is a descriptive sandbox. Time points can be serially dependent and the embeddings are not biologically validated. The map is only a summary of selected measurements/distances, not a statistical inference result.

## Methodology and interpretation

- **Generator:** Each feature is `baseline + signed_amplitude × (sigmoid(time) − 0.5) + Gaussian noise`, with resulting values clipped to `[0,1]`. Midpoint, baseline, sign, amplitude, width, and noise seed are recorded per feature. Gaussian noise is simulated independently per feature. Simulator time is an illustrative 0–5000 scale, not an observed patient time axis.
- **PCA:** Rows represent sampled time points; columns represent features. Each non-constant feature is centered and scaled *within the selected window*. A constant feature is removed from PCA, not imputed or hidden from correlation. PC1/PC2 scores refer to **time points**; loadings refer to **features**. PCA axes can flip sign without changing the mathematics.
- **Correlation:** Correlations are computed **between features across time points** within the same window. Spearman is Pearson correlation of average-tie ranks. When enabled, Gaussian weighting uses `w(t) = exp(-0.5 × ((t − focus)/bandwidth)^2)` and normalized weighted covariances. The *effective sample size* displayed is `(sum(w)^2 / sum(w^2))`, not a p-value or an independence-adjusted sample size. Constant features return `NA` correlation. All correlation matrices are symmetric and bounded by `[-1,1]`.
- **Interpretation:** Neighboring time points may be correlated, generated signals are designed to share structure, and clipping and standardization change estimated relationships. Correlations and PCA are **descriptive** here and should not be presented as statistical significance or biological validation.
- **Reproducibility:** Fixing the seed and generator parameters reproduces the simulator output under a compatible R random-number implementation. Download the settings and per-feature metadata; for rigorously reproducible deployments, also record the R version and lock installed package versions.

See [methodology notes](docs/METHODOLOGY.md) for technical details and [quick Korean UI smoke test](docs/SMOKE_TEST_KO.md) before deploying.

## GitHub Desktop guide (한국어)

1. 이 프로젝트 ZIP을 압축 해제하고 `trajectory-lab` 폴더를 연다.
2. R 또는 RStudio에서 `source("install.R")` 후 `shiny::runApp(".")`로 앱을 확인한다.
3. **GitHub Desktop → File → Add local repository**로 압축을 푼 폴더를 선택한다. Git 저장소가 아니라는 안내가 나오면 **Create a repository here**를 선택한다.
4. 첫 Commit 메시지는 `Initial Trajectory Lab release` 정도로 적고 **Publish repository**를 클릭한다.
5. GitHub에서 `README.md`가 자동 표시된다. 기능을 바꾸면 Commit과 Push를 통해 수정 이력을 남긴다.
6. 웹에 공개해서 실제로 실행시키려면 [DEPLOY.md](DEPLOY.md)를 참고한다. **GitHub에 소스를 올리는 것과 Shiny 앱을 호스팅하는 것은 별개**다.

> GitHub Desktop에서 기존 폴더를 대상으로 New repository를 만들 때 이미 `.git` 폴더가 있다면 중복 생성하지 말고 Add local repository를 이용한다. 이 ZIP에는 `.git`이 포함되지 않는다.

## Deployment

Use [Posit Connect Cloud](https://connect.posit.cloud/) for GitHub-connected deployment. A hosting-specific `manifest.json` must first be generated from a working R environment with `rsconnect::writeManifest()`. Full steps: [DEPLOY.md](DEPLOY.md).

No live instance has been deployed or verified as part of this release.

## Project layout

```text
trajectory-lab/
├── app.R                  # Shiny UI and server
├── R/
│   ├── experiments.R     # progressive PCA + mean/variance transformations
│   ├── plots_experiments.R # additional plots
│   ├── ordination.R      # PCA / distance matrix / PCoA / NMDS / Sammon
│   ├── plots_ordination.R # four-way dashboard, diagnostics, distance heatmap
│   ├── simulation.R       # deterministic trajectory generator
│   ├── analysis.R         # PCA and weighted/unweighted correlation
│   ├── io.R               # validated CSV import and export notes
│   └── plots.R            # chart rendering functions
├── tests/run_tests.R      # dependency-free computational checks
├── .github/workflows/     # test workflow
├── www/                   # responsive stylesheet + icon
├── docs/METHODOLOGY.md
├── docs/SMOKE_TEST_KO.md
├── examples/demo-beta.csv  # synthetic upload example
├── install.R
├── DEPLOY.md
├── CHANGELOG.md
├── VERSION
└── LICENSE
```

## License

[MIT License](LICENSE). © 2026 ka27cari. Includes **no original proprietary research data, RDS objects, or user-specific paths**.
