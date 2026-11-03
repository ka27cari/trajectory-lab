# Methodology and limitations

This document is not a scientific methods section for any empirical study. It describes the math used by an interactive demonstration.

## Synthetic beta-like time courses

For time `t` and feature `j`:

```
p_j(t) = 1 / (1 + exp(-4 × (t - midpoint_j) / width_j))
y_j(t) = clamp(baseline_j + amplitude_j × (p_j(t) - 0.5)
               + Normal(0, noise_sd_j), 0, 1)
```

Negative amplitudes produce downward shifts, positive amplitudes upward shifts. The simulator randomly varies parameters within the user-specified preset and records each feature's exact generated parameters in metadata. The noise seed is deterministic and independent for each feature. This is a visual model of a bounded signal, **not** a calibrated methylation assay noise model; clipping generates a point mass at the bounds in some settings.

## PCA definition

For the selected observation window, let `X` be `N time points × P features`. Column means and standard deviations are estimated from **the selected N rows**, not a broader dataset. Constant columns are excluded. `stats::prcomp(X, center=TRUE, scale.=TRUE)` computes singular vectors, time-point score coordinates and variable loadings. Eigenvalue fractions are `(sdev^2) / sum(sdev^2)`.

The displayed path joins score points in ascending time order. It does not represent causal evolution or a fitted dynamical model. PC signs are arbitrary; loadings and scores must be interpreted consistently.

## Progressive PCA and feature loading walks

PCA is recomputed from all active, non-constant features in the selected window to define a **reference** 2D coordinate system. When adding features in *fixed-basis* mode, the partial score for time `i` is `sum_{j in included} z_ij × loading_j,1:2`, where each `z_ij` uses the full-window means and SDs and the loadings come from full-window PCA. The final partial score equals the full PC1/PC2 score.

In *refit* mode, PCA is recomputed for each new feature set. A two-dimensional orthogonal Procrustes transform aligns the step's scores with the full reference by minimizing squared Euclidean differences. This **only improves visual orientation and does not make refitted principal directions equal or remove information loss**. Results should not be read as tracking an invariant PC coordinate system.

The loading walk accumulates full-reference PC1/PC2 loading vectors in the chosen feature order. This is a visual construction, not a score trajectory or proof of causality. For uploaded CSVs, missing transition midpoint metadata fall back to creation order.

## Means, variance and controlled transformations

The active sample window is split into two consecutive halves by row order, not by two independent biological conditions. Mean change and late/early variance ratio are descriptive. The variance ratio is undefined when early variance is essentially zero. A log2 ratio of −1 in the scatter means half as much variance in the late half; a log2 ratio of 0 means equal variances.

The experiment tab applies `y = clamp(mean(x) + offset + multiplier * (x − mean(x)), 0, 1)` to exactly one feature, in the selected time window only. The transform does not add random noise and does not mutate the baseline dataset. When the apply toggle is on, PCA, correlations and active-window exports use the altered values; generator metadata and full-data exports remain unchanged. Clamping can change the requested mean shift and variance multiplication. The UI shows realized values instead of pretending that requested targets were always met.

## Distance-aware ordination (v1.2.0)

For time-point ordination, the observation matrix has time points on rows and non-constant beta-like features on columns. For feature ordination, the same preprocessed matrix is transposed, so each point corresponds to one full feature profile across time. Standardization is optional (the default z-scores each feature within the analysis window). For display and computation time, up to 200 ordered time points are evenly sampled **only in the distance workbench**; sampling occurs after selecting the observation window.

Supported dissimilarities are Euclidean (L2), Manhattan (L1), Canberra (sum of absolute coordinate differences divided by sums of coordinate magnitudes, with R `dist` zero conventions), Bray–Curtis (sum absolute differences divided by sum of paired nonnegative values; unavailable for negative standardized inputs or all-zero pairs), and correlation dissimilarity (1 − Pearson r computed across the columns of each pair of object profiles). Correlation distance requires nonconstant profiles and is not guaranteed to be Euclidean.

PCA is calculated on the preprocessed measurement matrix, irrespective of the distance selection. PCoA uses `stats::cmdscale(distance, eig = TRUE, k = 2)` with no additive correction and therefore may have negative eigenvalues for non-Euclidean dissimilarities. NMDS uses `MASS::isoMDS` and Sammon mapping uses `MASS::sammon`, initialized by classical MDS. Their stress values optimize distinct loss functions and are **not directly comparable** to PCA variance explained.

The common rank-diagnostic is `cor(original pair dissimilarities, final 2D Euclidean distances, method = 'spearman')`. The common relative distance error is `sqrt(sum((d - alpha*r)^2) / sum(d^2))`, where `r` is the pair distance in the fitted 2D map and `alpha = sum(d*r)/sum(r^2)` (0 when denominator is numerically zero). It is a scale-adjusted summary, **not** an algorithm-specific optimized loss. Missing or degenerate methods are reported independently; zero off-diagonal distances prohibit NMDS/Sammon instead of being silently nudged with random noise.

PCA/PCoA geometries can be equivalent up to rigid rotation/reflection when the selected distance is Euclidean and both use the same centered measurement matrix. This property is tested directly using interpoint distances rather than coordinate equality, which would be misleading because axis orientation is arbitrary.

## Correlation definition

For feature pair `x` and `y`, optional nonnegative Gaussian time weights `w_i`, and their weighted means `mu_x`, `mu_y`:

```
cov_w(x,y) = sum_i(w_i × (x_i - mu_x) × (y_i - mu_y)) / sum_i(w_i)
r_w(x,y) = cov_w(x,y) / sqrt(cov_w(x,x) × cov_w(y,y))
```

When Gaussian weighting is off, all weights equal one. Spearman correlation replaces each variable with its average-tie ranks before applying the same formula. Correlation for constant variables is `NA`, not zero. The `effective_n = (sum w)^2 / (sum w^2)` statistic reflects concentration of weights, **not independent biological sample count**.

**Change from the legacy prototype:** Correlation is no longer indexed by time when its axes describe feature pairs; variance terms are correctly derived from their respective variables; a single symmetric, normalized weight vector is used per observation; a valid coefficient is bounded in `[-1,1]`. The old ad hoc balance / moving-average / adjacency-power switches are intentionally not reproduced because the previous implementations mixed concepts or lacked a defensible statistical interpretation.

## Data model

All features must be finite numeric beta-scale signals between 0 and 1. `time` must be finite, strictly increasing, and unique. The upload pathway rejects missing values instead of silently imputing them. The UI caps CSV sizes and matrix dimensions. Generated and uploaded matrices follow the same analysis path.

## Limitations

- Synthetic time points do not represent independent subjects.
- Pearson and Spearman coefficients are descriptive, and the app has no multiple-comparison testing, confidence intervals, or significance thresholds.
- No normalization for EPIC assay variants or probe filtering is implemented.
- No batch effect, clinical covariate adjustment, deconvolution, missing-value imputation, or longitudinal mixed modeling is implemented.
- Input CSV may contain real measurements, but importing data does not make the analysis validated for any specific study.
- Browser visualizations and cloud service availability require an actual deployed Shiny instance; GitHub source code alone is not a hosted service.
