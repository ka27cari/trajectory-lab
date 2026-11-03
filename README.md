# Trajectory Lab 🧬

**An interactive R Shiny workspace for exploring synthetic methylation-like trajectories, PCA, distance-based embeddings, and the difference between mean shifts and changes in spread.**

**R Shiny** · **PCA / PCoA / NMDS / Sammon** · **5 distance definitions** · **CSV import/export**

Trajectory Lab makes it easier to ask *what changed, and why does the map look different?* Adjust simulated time courses, change the set of features entering PCA, modify a feature's mean or spread, and compare different ways of representing relationships among time points or features.

Trajectory Lab is an exploratory research and teaching application built around simulated data. It does not include patient-level data.

## At a glance

| Experiment | Explore | Interpret |
|:--|:--|:--|
| **Trajectories** | Generate bounded synthetic signals or upload a beta-scale CSV | Changes across the observation window |
| **Progressive PCA** | Add features step by step, refit/align the map or accumulate fixed-axis contributions | Why PC1/PC2 point positions change as inputs change |
| **Mean & variance** | Shift a feature's mean or compress/expand its spread | Whether location, variability or both changed |
| **Distance maps** | Compare PCA, PCoA, NMDS and Sammon in one 2×2 view | How the distance definition affects geometry |
| **Correlations** | Inspect Pearson/Spearman matrices and individual feature pairs | Similarity of feature trajectories |

### One dataset, four views

The **Distance maps** tab lays out **PCA, PCoA, NMDS and Sammon** side by side. Try Euclidean, Manhattan, Canberra, Bray–Curtis or `1 − Pearson correlation` to see what changes in the distance-based methods. Switch between mapping **time points** and mapping **feature profiles**.

**Important:** PCA works on the feature-value matrix; it **does not** use the distance selector. PCoA, NMDS and Sammon use the selected pairwise dissimilarities. Axes are not directly comparable between methods, so inspect geometric relationships and the accompanying diagnostics rather than absolute coordinates. [Read the methodology](docs/METHODOLOGY.md).

## Run locally

Requirements: R, with `shiny`, `bslib`, `ggplot2` and `MASS` packages. From the repository's **top-level folder** (the folder that contains `app.R`):

```r
source("install.R")
shiny::runApp(".")
```

In a terminal with R installed, computational checks can be run with:

```sh
Rscript tests/run_tests.R
```

Automated computational checks are configured in `.github/workflows/`; run the app locally to evaluate interactive behavior.

## Try these experiments

1. **Follow PC1/PC2 as features are added.** Under **PCA → Feature-by-feature PCA**, step through the input variables. Choose between a fresh PCA fit (aligned for display) and cumulative contributions to a fixed PCA basis. The loading-vector walk is *not* a time-point score path.
2. **Move the mean without deliberately changing spread.** Under **Mean & variance**, set a nonzero mean offset and a spread multiplier of `1`. Then set the offset to `0` and multiplier to `0.5` to narrow variability. Boundary clipping to `[0, 1]` may affect both quantities; check the *realized* statistics.
3. **Change what “close” means.** In **Distance maps**, start with Euclidean distance, then select Manhattan or correlation distance and click **Recalculate four maps**. Compare the overview and quality diagnostics. For Bray–Curtis, disable z-score standardization because signed inputs are inappropriate.
4. **Inspect a correlation pair.** Under **Correlations**, click a matrix cell to inspect the two time courses, or use Pearson/Spearman and optional time weighting.

## Bring your own data

Upload a CSV with one numeric `time` column and **2–40** numeric feature columns containing beta-scale values between **0 and 1**. Time values must be strictly increasing, with no missing or nonfinite cells. The import accepts **5–1000 rows** and a **3 MB** maximum file size.

```csv
time,CpG_01,CpG_02,CpG_03
0,0.42,0.73,0.29
100,0.44,0.71,0.31
200,0.46,0.69,0.32
300,0.48,0.65,0.34
400,0.50,0.62,0.36
```

A full synthetic example is in [`examples/demo-beta.csv`](examples/demo-beta.csv). User uploads are processed in the active Shiny session by this application's code, which does not explicitly save them to disk or a database; hosting infrastructure may have its own policies. Do **not** upload sensitive research data to an unverified public deployment.

## Results and reproducibility

Downloadable outputs include input/selected-window matrices, simulation metadata, PCA scores, progressive PCA results, correlations, distance matrices, embedding coordinates, map quality diagnostics, settings and PNG plots. Fix the random seed to reproduce the built-in simulation in compatible R environments.

Method-specific notes:

- **Progressive PCA:** Refit-and-align mode recomputes PCA with a growing set of features; fixed-basis mode accumulates contributions using a shared reference basis. Neither means that “adding PC3” changes already fitted PC1/PC2 scores.
- **Mean/spread experiments:** Time-window halves are consecutive observations, **not independent patient groups**. Clipping may alter the requested transformation.
- **Distance preservation:** PCoA, NMDS and Sammon optimize different objectives. Correlation of distance ranks, scaled distance error and each method's stress must be interpreted separately.
- **Scientific scope:** Illustrative simulation and descriptive ordination, **not** EWAS, deconvolution, batch correction, statistical significance, causal inference or clinical prediction.

For definitions, caveats and model assumptions, see [Methodology](docs/METHODOLOGY.md).

## Project structure

```text
trajectory-lab/
├── app.R                 # Shiny UI and server
├── R/                    # generator, analyses, ordinations, plotting, IO
├── www/                  # app stylesheet and icon
├── examples/             # synthetic demonstration CSV
├── tests/                # R computational tests
├── docs/                 # methodology, smoke test, overview illustration
├── .github/workflows/    # CI
├── install.R
├── DEPLOY.md
├── CHANGELOG.md
└── README.md
```

## Online demo and deployment

This repository provides the Shiny source code. Running the app online requires an R Shiny hosting service. See [Deployment](DEPLOY.md) for the hosting workflow. No public demo URL is currently provided.

## License

[MIT](LICENSE) · © 2026 ka27cari.
