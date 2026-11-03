# Changelog

## 1.2.0 | Distance-aware 2D embedding workbench

- Continued the interrupted v1.1.0 source instead of restarting from v1.0.1; retained feature-step PCA, controlled mean/spread experiments and the correlation-matrix hotfix.
- Added a one-screen comparison of standardized or unscaled PCA, classical PCoA, nonmetric MDS and Sammon mapping, with Euclidean, Manhattan, Canberra, Bray-Curtis, and 1 − Pearson correlation dissimilarities.
- Added time-point/feature orientation, deterministic time-point subsampling (up to 200), pairwise distance heatmap, relative distance error and pair-distance rank correlation, an individual map and a diagnostics scatter.
- Added downloadable all-method coordinates, quality diagnostics, distance matrix and overview PNG; expanded source parsing and numerical tests, including PCA–PCoA Euclidean equivalence.
- Restored one-feature starting point for fixed-basis PCA, while refitting correctly starts with two. Reduced chart legend widths, tick counts and long numeric tick labels to minimize collisions.
- PCoA, NMDS and Sammon validate each method separately and show readable unavailable messages for degenerate inputs. MASS is the only new dependency.
- **Execution note:** R and browser testing are not available in this environment. Run the provided tests and smoke-test checklist before claiming a validated deployment.

## 1.1.0 | PCA step-through and mean/variance experiments

- Brought back stepwise visual exploration from the legacy `dimReduction.R` idea: growing-feature PCA refits (orthogonally aligned for visual comparison) and an explicitly different fixed-basis cumulative loading walk.
- Added two modes for sequential PCA and three feature inclusion orders.
- Added early/late mean–variance comparisons across all features, plus a controlled mean-offset/spread-multiplier experiment.
- Experimental transformation optionally propagates to PCA, correlations, and active-window exports, while preserving original raw full-data downloads.
- Added separate reports for active-window, mean–variance and progressive PCA results; progressive PCA figures can be downloaded as PNG.
- Reworked color-bar sizes and tick breaks; reduced heatmap axis labels for high feature counts, leaving all matrix cells available.
- Added deterministic computational regression tests for refit/fixed PCA, cumulative loadings, mean shifts, variance compression and clipping. CI now attempts Shiny app object construction after dependency installation.
- **Execution note:** tests are written but could not be run locally because no R interpreter is available; the GitHub Actions workflow runs them after a push.

## 1.0.1 | Correlation matrix hotfix

- Preserve the correlation matrix dimensions when bounding computed coefficients to [-1, 1]. The previous scalar-first `pmin()`/`pmax()` calls could drop matrix attributes and cause the `dimnames()` error.
- Restore the correlation heatmap, correlation CSV export, pair statistics, and `MEAN |r|` metric that depended on the failed calculation.
- Add regression cases covering Pearson and Spearman, unweighted and Gaussian-weighted paths, two-feature datasets, matrix dimensions, names, symmetry, and finite mean absolute correlation.
- Calculation test suite remains **not executed locally** because this environment has no R interpreter. GitHub Actions will execute it once pushed; this is a code-level fix, not an execution-verified release.

## 1.0.0 | Initial standalone service

- Rebuilt legacy time-course/PCA and correlation prototypes as one Shiny app.
- Added deterministic, configurable, bounded synthetic signals with four timing presets.
- Added custom feature creation/undo and strict beta-scale CSV importing.
- Added shared observation-window filtering, PCA, scree plots, loadings, Pearson/Spearman and optional Gaussian time weighting.
- Added clickable correlation heatmap, pair inspector and six data/figure downloads.
- Isolated computation from UI and introduced base-R validation tests and GitHub Actions CI.
- Only a small synthetic example CSV is bundled; no empirical research data, local paths, platform-specific dependencies or external research results are included.
