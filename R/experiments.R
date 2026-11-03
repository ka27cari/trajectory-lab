# Pure analysis functions for the PCA step-through and shift/spread lab.
# No Shiny, ggplot, or new dependencies; each function is testable in base R.

active_features <- function(data) {
  features <- feature_names(data)
  X <- as.matrix(data[, features, drop = FALSE])
  sds <- apply(X, 2L, stats::sd)
  features[is.finite(sds) & sds > 1e-10]
}

feature_sequence <- function(data, metadata = NULL, mode = "created") {
  eligible <- active_features(data)
  if (length(eligible) < 2L) {
    stop("At least two non-constant features are needed.", call. = FALSE)
  }
  if (identical(mode, "created")) return(eligible)
  if (identical(mode, "late")) {
    midpoints <- if (is.null(metadata)) rep(NA_real_, length(eligible)) else
      metadata$midpoint[match(eligible, metadata$feature)]
    if (all(is.na(midpoints))) return(eligible)
    return(eligible[order(is.na(midpoints),
                          -ifelse(is.na(midpoints), 0, midpoints),
                          seq_along(eligible))])
  }
  if (identical(mode, "loading")) {
    fit <- stats::prcomp(data[, eligible, drop = FALSE], center = TRUE, scale. = TRUE)
    weights <- abs(fit$rotation[, 1L])
    return(eligible[order(-weights[eligible], seq_along(eligible))])
  }
  stop("Unknown feature ordering.", call. = FALSE)
}

progressive_pca <- function(data, sequence, step, mode = "refit") {
  eligible <- active_features(data)
  if (length(eligible) < 2L || length(sequence) != length(eligible) ||
      anyDuplicated(sequence) || !setequal(sequence, eligible) ||
      length(step) != 1L || is.na(step) || step != as.integer(step) ||
      step < 1L || step > length(sequence) ||
      !mode %in% c("refit", "fixed")) {
    stop("Invalid progressive PCA settings.", call. = FALSE)
  }
  if (identical(mode, "refit") && step < 2L) {
    stop("Refitted PCA needs at least two included features.", call. = FALSE)
  }
  X <- as.matrix(data[, eligible, drop = FALSE])
  full_fit <- stats::prcomp(X, center = TRUE, scale. = TRUE)
  reference <- full_fit$x[, 1:2, drop = FALSE]
  selected <- sequence[seq_len(step)]
  if (identical(mode, "fixed")) {
    # Partial sums in ONE fixed reference basis; not a PCA refit.
    Z <- scale(X, center = full_fit$center, scale = full_fit$scale)
    scores <- Z[, selected, drop = FALSE] %*%
      full_fit$rotation[selected, 1:2, drop = FALSE]
  } else {
    # Refit PCA on a growing variable set. Orthogonal Procrustes alignment
    # removes arbitrary rotations/reflections for graphical comparison only.
    subset_fit <- stats::prcomp(X[, selected, drop = FALSE],
                                center = TRUE, scale. = TRUE)
    unaligned <- subset_fit$x[, 1:2, drop = FALSE]
    sv <- svd(crossprod(unaligned, reference))
    scores <- unaligned %*% (sv$u %*% t(sv$v))
  }
  colnames(scores) <- c("PC1", "PC2")
  colnames(reference) <- c("PC1", "PC2")
  list(current = data.frame(time = data$time, scores, check.names = FALSE),
       reference = data.frame(time = data$time, reference, check.names = FALSE),
       sequence = sequence, included = selected, step = step, mode = mode)
}

loading_walk <- function(data, sequence, step) {
  eligible <- active_features(data)
  if (length(sequence) != length(eligible) || anyDuplicated(sequence) ||
      !setequal(sequence, eligible) || length(step) != 1L || is.na(step) ||
      step != as.integer(step) || step < 1L || step > length(sequence)) {
    stop("Invalid loading walk settings.", call. = FALSE)
  }
  fit <- stats::prcomp(data[, eligible, drop = FALSE], center = TRUE, scale. = TRUE)
  shifts <- fit$rotation[sequence[seq_len(step)], 1:2, drop = FALSE]
  coords <- rbind(c(0, 0), apply(shifts, 2L, cumsum))
  colnames(coords) <- c("PC1", "PC2")
  data.frame(step = 0:step,
             feature = c("Origin", sequence[seq_len(step)]),
             PC1 = coords[, 1L], PC2 = coords[, 2L],
             stringsAsFactors = FALSE)
}

half_window_summary <- function(data) {
  if (nrow(data) < 5L) stop("At least five time points are required.", call. = FALSE)
  early <- seq_len(floor(nrow(data) / 2))
  late <- seq.int(max(early) + 1L, nrow(data))
  features <- feature_names(data)
  make_row <- function(feature) {
    a <- data[[feature]][early]
    b <- data[[feature]][late]
    va <- stats::var(a)
    vb <- stats::var(b)
    data.frame(feature = feature, mean_early = mean(a), mean_late = mean(b),
               delta_mean = mean(b) - mean(a), var_early = va, var_late = vb,
               variance_ratio = if (va < 1e-12) NA_real_ else vb / va,
               stringsAsFactors = FALSE)
  }
  do.call(rbind, lapply(features, make_row))
}

perturb_feature <- function(data, feature, shift = 0, spread = 1) {
  if (length(feature) != 1L || is.na(feature) ||
      !feature %in% feature_names(data) ||
      length(shift) != 1L || length(spread) != 1L ||
      !is.finite(shift) || !is.finite(spread) ||
      shift < -1 || shift > 1 || spread < 0 || spread > 4) {
    stop("Invalid shift/spread experiment settings.", call. = FALSE)
  }
  baseline <- data[[feature]]
  proposed <- mean(baseline) + shift + spread * (baseline - mean(baseline))
  altered <- pmax(0, pmin(1, proposed))
  transformed <- data
  transformed[[feature]] <- altered
  list(data = transformed,
       paired = data.frame(time = data$time, baseline = baseline, altered = altered),
       actual_shift = mean(altered) - mean(baseline),
       variance_ratio = if (stats::var(baseline) < 1e-12) NA_real_ else
         stats::var(altered) / stats::var(baseline),
       clipped = sum(proposed < 0 | proposed > 1),
       total = length(proposed))
}
