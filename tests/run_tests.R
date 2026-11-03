# Base R tests: no external testing framework is required.
# From repository root: Rscript tests/run_tests.R

source("R/simulation.R")
source("R/analysis.R")
source("R/experiments.R")
source("R/ordination.R")
source("R/io.R")

expect_error <- function(expression) {
  result <- try(force(expression), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}

# Everything must parse (including UI files, without executing UI dependencies).
for (file in c("app.R", "R/simulation.R", "R/analysis.R",
               "R/io.R", "R/plots.R", "R/experiments.R",
               "R/plots_experiments.R", "R/ordination.R", "R/plots_ordination.R")) {
  parse(file = file)
}

# Same seed and settings always yield identical data.
a <- simulate_dataset(n_features = 8, n_timepoints = 90, seed = 123)
b <- simulate_dataset(n_features = 8, n_timepoints = 90, seed = 123)
stopifnot(identical(a, b), nrow(a$data) == 90,
          length(feature_names(a$data)) == 8,
          all(as.matrix(a$data[, -1L]) >= 0),
          all(as.matrix(a$data[, -1L]) <= 1))

# The simulator must not mutate the caller's global random stream.
set.seed(888)
old_rng <- .Random.seed
invisible(simulate_dataset(seed = 24))
stopifnot(identical(old_rng, .Random.seed))

expect_error(simulate_dataset(n_features = 1))
expect_error(simulate_dataset(noise_sd = -1))
expect_error(simulate_dataset(pattern = "unknown"))

# Custom feature can be added and undone without changing original profiles.
with_custom <- add_custom_feature(a, 1500, -0.5, 0.55, 700, 0.04)
stopifnot(ncol(with_custom$data) == ncol(a$data) + 1L,
          nrow(with_custom$metadata) == nrow(a$metadata) + 1L,
          identical(undo_custom_feature(with_custom), a))
expect_error(undo_custom_feature(a))

# PCA uses time rows as observations and features as columns.
pca <- run_pca(a$data)
stopifnot(nrow(pca$scores) == nrow(a$data),
          nrow(pca$loadings) == length(feature_names(a$data)),
          abs(sum(pca$explained$fraction) - 1) < 1e-10,
          all(is.finite(as.matrix(pca$scores[, -1L]))))

# Constant features are excluded from PCA without destroying correlation output.
constant <- a$data
constant$constant <- 0.4
pca_constant <- run_pca(constant)
stopifnot("constant" %in% pca_constant$dropped)
expect_error(run_pca(data.frame(time = 1:9, A = rep(0.5, 9), B = rep(0.5, 9))))

# Conventional Pearson/Spearman must match stats::cor with equal weights.
X <- as.matrix(a$data[, -1L])
reference_pearson <- stats::cor(X, method = "pearson")
reference_spearman <- stats::cor(X, method = "spearman")
test_pearson <- weighted_correlation(X, rep(1, nrow(X)), "pearson")
test_spearman <- weighted_correlation(X, rep(1, nrow(X)), "spearman")
stopifnot(isTRUE(all.equal(unname(reference_pearson), unname(test_pearson), tolerance = 1e-10)),
          isTRUE(all.equal(unname(reference_spearman), unname(test_spearman), tolerance = 1e-10)))

# Regression: correlation must remain a named two-dimensional matrix.
# v1.0.0 lost dimensions when clipping with pmin(1, pmax(-1, matrix)),
# causing a dimnames() error in the correlation tab and MEAN |r| = N/A.
for (method in c("pearson", "spearman")) {
  for (weighted_flag in c(FALSE, TRUE)) {
    out <- run_correlation(a$data, method = method, weighted = weighted_flag,
                           focus = 1000, bandwidth = 400)
    stopifnot(is.matrix(out$matrix),
              identical(dim(out$matrix), c(8L, 8L)),
              identical(rownames(out$matrix), feature_names(a$data)),
              identical(colnames(out$matrix), feature_names(a$data)),
              is.finite(out$mean_abs),
              out$mean_abs >= 0, out$mean_abs <= 1,
              isTRUE(all.equal(out$matrix, t(out$matrix), tolerance = 1e-12)))
  }
}
# Smallest valid feature count must also retain a 2x2 matrix.
two_features <- simulate_dataset(n_features = 2, n_timepoints = 20, seed = 19)
for (method in c("pearson", "spearman")) {
  two_cor <- run_correlation(two_features$data, method = method)
  stopifnot(is.matrix(two_cor$matrix),
            identical(dim(two_cor$matrix), c(2L, 2L)),
            identical(rownames(two_cor$matrix), feature_names(two_features$data)),
            is.finite(two_cor$mean_abs))
}

# Gaussian-weighted correlation should remain symmetric and in [-1,1].
weighted <- run_correlation(a$data, weighted = TRUE, focus = 900, bandwidth = 300)
stopifnot(isTRUE(all.equal(weighted$matrix, t(weighted$matrix))),
          max(abs(weighted$matrix), na.rm = TRUE) <= 1 + 1e-12,
          weighted$effective_n >= 1,
          weighted$effective_n <= nrow(a$data))

perfect <- cbind(A = 1:10, B = 11:20, C = 10:1, D = rep(2, 10))
perfect_cor <- weighted_correlation(perfect, rep(1, 10))
stopifnot(abs(perfect_cor["A", "B"] - 1) < 1e-12,
          abs(perfect_cor["A", "C"] + 1) < 1e-12,
          is.na(perfect_cor["A", "D"]))

# Short/empty observation windows should give clear errors.
expect_error(window_data(a$data, c(1, 1)))
expect_error(window_data(a$data, c(1, 2)))
stopifnot(nrow(window_data(a$data, c(0, 1000))) >= 5L)

# CSV input: round-trip and malformed inputs.
file <- tempfile(fileext = ".csv")
utils::write.csv(a$data, file, row.names = FALSE)
loaded <- read_beta_csv(file)
stopifnot(isTRUE(all.equal(loaded$data, a$data, tolerance = 1e-10)),
          identical(loaded$source, "uploaded"))

bad <- a$data
bad$CpG_01[1] <- 1.2
utils::write.csv(bad, file, row.names = FALSE)
expect_error(read_beta_csv(file))

bad <- a$data
bad$time[2L] <- bad$time[1L]
utils::write.csv(bad, file, row.names = FALSE)
expect_error(read_beta_csv(file))

bad <- a$data
bad$CpG_01[1L] <- NA_real_
utils::write.csv(bad, file, row.names = FALSE)
expect_error(read_beta_csv(file))

# Feature-by-feature PCA: fixed-basis partial sums must reproduce full scores.
sequence <- feature_sequence(a$data, a$metadata, "created")
stopifnot(identical(sequence, feature_names(a$data)),
          identical(sort(feature_sequence(a$data, a$metadata, "late")), sort(sequence)),
          identical(sort(feature_sequence(a$data, a$metadata, "loading")), sort(sequence)))
fixed_full <- progressive_pca(a$data, sequence, length(sequence), mode = "fixed")
stopifnot(isTRUE(all.equal(fixed_full$current, fixed_full$reference,
                           tolerance = 1e-9)))
fixed_one <- progressive_pca(a$data, sequence, 1L, mode = "fixed")
stopifnot(nrow(fixed_one$current) == nrow(a$data),
          all(is.finite(as.matrix(fixed_one$current[, -1L]))))
refit_two <- progressive_pca(a$data, sequence, 2L, mode = "refit")
stopifnot(all(is.finite(as.matrix(refit_two$current[, -1L]))),
          identical(refit_two$included, sequence[1:2]))
refit_full <- progressive_pca(a$data, sequence, length(sequence), mode = "refit")
stopifnot(isTRUE(all.equal(refit_full$current, refit_full$reference,
                           tolerance = 1e-9)))
walk <- loading_walk(a$data, sequence, 4L)
stopifnot(nrow(walk) == 5L, identical(walk$step, 0:4),
          all(is.finite(as.matrix(walk[, c("PC1", "PC2")]))))
stopifnot(length(active_features(constant)) == 8L,
          length(feature_sequence(constant, NULL, "created")) == 8L)
expect_error(progressive_pca(a$data, sequence, 1L, mode = "refit"))
expect_error(progressive_pca(a$data, sequence, 100L, mode = "fixed"))
expect_error(feature_sequence(a$data, a$metadata, "invalid"))

# Means and spread: separate location changes from dispersion changes.
controlled <- data.frame(time = 1:8,
                         A = c(0.2, 0.3, 0.4, 0.5, 0.4, 0.45, 0.5, 0.55),
                         B = c(0.5, 0.5, 0.5, 0.5, 0.6, 0.6, 0.6, 0.6))
shifted <- perturb_feature(controlled, "A", shift = 0.1, spread = 1)
compressed <- perturb_feature(controlled, "A", shift = 0, spread = 0.5)
stopifnot(abs(shifted$actual_shift - 0.1) < 1e-10,
          abs(shifted$variance_ratio - 1) < 1e-10,
          abs(compressed$actual_shift) < 1e-10,
          abs(compressed$variance_ratio - 0.25) < 1e-10,
          isTRUE(all.equal(controlled, perturb_feature(controlled, "A")$data,
                           tolerance = 1e-12)),
          identical(controlled$B, shifted$data$B))
clipped <- perturb_feature(controlled, "A", shift = 0.9, spread = 2)
stopifnot(all(clipped$data$A >= 0 & clipped$data$A <= 1),
          clipped$clipped > 0L)
expect_error(perturb_feature(controlled, "missing"))
expect_error(perturb_feature(controlled, "A", spread = -1))
changes <- half_window_summary(controlled)
stopifnot(nrow(changes) == 2L,
          abs(changes$delta_mean[changes$feature == "A"] -
              (mean(tail(controlled$A, 4L)) - mean(head(controlled$A, 4L)))) < 1e-12,
          is.na(changes$variance_ratio[changes$feature == "B"]))

# Distance explorer: inputs, geometry, dataset orientation and valid metrics.
prep <- prepare_ordination(a$data, orientation = "time", standardize = TRUE,
                           max_points = 20L)
stopifnot(prep$n_objects == 20L, prep$total_objects == 90L,
          length(prep$time) == 20L, prep$time[1L] == a$data$time[1L],
          tail(prep$time, 1L) == tail(a$data$time, 1L),
          nrow(prep$X) == 20L, ncol(prep$X) == 8L)
for (metric in c("euclidean", "manhattan", "canberra", "correlation")) {
  D <- make_dissimilarity(prep, metric)
  mat <- as.matrix(D)
  stopifnot(inherits(D, "dist"), nrow(mat) == 20L,
            isTRUE(all.equal(mat, t(mat), tolerance = 1e-10)),
            all(diag(mat) == 0),
            all(is.finite(mat)), all(mat >= 0))
}
stopifnot(isTRUE(all.equal(unname(as.numeric(make_dissimilarity(prep, "manhattan"))),
                           unname(as.numeric(stats::dist(prep$X, method = "manhattan"))),
                           tolerance = 1e-12)))
raw_prep <- prepare_ordination(a$data, standardize = FALSE, max_points = 20L)
bray <- as.matrix(make_dissimilarity(raw_prep, "bray"))
stopifnot(all(bray >= 0 & bray <= 1),
          abs(bray[1, 2] - sum(abs(raw_prep$X[1, ] - raw_prep$X[2, ])) /
                sum(raw_prep$X[1, ] + raw_prep$X[2, ])) < 1e-12)
expect_error(make_dissimilarity(prep, "bray"))
expect_error(make_dissimilarity(prep, "not-a-distance"))
features_prep <- prepare_ordination(a$data, orientation = "feature",
                                   standardize = TRUE, max_points = 20L)
stopifnot(features_prep$n_objects == 8L,
          features_prep$total_objects == 8L,
          ncol(features_prep$X) == 90L,
          identical(features_prep$labels, feature_names(a$data)))
expect_error(prepare_ordination(a$data, max_points = 1L))

# On a standardized Euclidean distance matrix PCA and classical PCoA must
# produce equivalent interpoint distances in the first two dimensions.
suite <- run_ordination_suite(a$data, distance = "euclidean",
                              standardize = TRUE, max_points = 20L)
stopifnot(identical(names(suite$fits), ordination_methods),
          nrow(suite$summary) == 4L,
          suite$summary$status[1L] == "OK",
          suite$summary$status[2L] == "OK",
          nrow(ordination_scores(suite)) >= 40L,
          nrow(ordination_pair_data(suite, "PCA")) == 190L,
          identical(dim(ordination_distance_preview(suite, 9L)), c(9L, 9L)))
pca_xy <- suite$fits$PCA$scores[, c("Dim1", "Dim2")]
pcoa_xy <- suite$fits$PCoA$scores[, c("Dim1", "Dim2")]
stopifnot(isTRUE(all.equal(as.numeric(stats::dist(pca_xy)),
                           as.numeric(stats::dist(pcoa_xy)), tolerance = 1e-6)))
for (method in ordination_methods) {
  obj <- suite$fits[[method]]
  if (is.null(obj$error)) {
    stopifnot(nrow(obj$scores) == 20L,
              all(is.finite(as.matrix(obj$scores[, c("Dim1", "Dim2")]))),
              obj$quality$distance_error >= 0,
              obj$quality$distance_error <= 1 + 1e-10,
              is.na(obj$quality$distance_rank) ||
                abs(obj$quality$distance_rank) <= 1 + 1e-10)
  } else {
    stopifnot(is.character(obj$error), nzchar(obj$error))
  }
}

# Duplicates in the mapped objects are a controlled limitation for MASS methods.
with_duplicates <- a$data
with_duplicates[2L, -1L] <- with_duplicates[1L, -1L]
with_duplicates[3L, -1L] <- with_duplicates[1L, -1L]
duplicate_maps <- run_ordination_suite(with_duplicates, standardize = FALSE,
                                       max_points = 90L)
stopifnot(!is.null(duplicate_maps$fits$NMDS$error),
          !is.null(duplicate_maps$fits$Sammon$error),
          is.null(duplicate_maps$fits$PCA$error))

# No accidental mutation of the original dataset by distance operations.
stopifnot(identical(a, b))
unlink(file)
cat("PASS: simulator, reproducibility, PCA, correlations, CSV validation, progressive PCA, mean/variance, pairwise distances, PCA–PCoA equivalence, method fallbacks, and source parsing.\n")
