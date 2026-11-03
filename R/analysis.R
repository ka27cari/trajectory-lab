feature_names <- function(data) {
  setdiff(names(data), "time")
}

window_data <- function(data, window) {
  if (length(window) != 2L || any(!is.finite(window)) ||
      window[1] >= window[2]) {
    stop("Choose a valid observation window.", call. = FALSE)
  }
  selected <- data[data$time >= window[1] & data$time <= window[2], , drop = FALSE]
  if (nrow(selected) < 5L) {
    stop("The selected window needs at least five time points.", call. = FALSE)
  }
  selected
}

run_pca <- function(data) {
  features <- feature_names(data)
  X <- as.matrix(data[, features, drop = FALSE])
  sds <- apply(X, 2L, stats::sd)
  keep <- is.finite(sds) & sds > 1e-10
  if (sum(keep) < 2L) {
    stop("PCA requires at least two non-constant features.", call. = FALSE)
  }

  fit <- stats::prcomp(X[, keep, drop = FALSE], center = TRUE, scale. = TRUE)
  if (ncol(fit$x) < 2L) {
    stop("PCA requires at least two components.", call. = FALSE)
  }
  variance <- fit$sdev^2
  fraction <- variance / sum(variance)
  scores <- data.frame(time = data$time, PC1 = fit$x[, 1L], PC2 = fit$x[, 2L])
  loadings <- data.frame(feature = rownames(fit$rotation),
                         PC1 = fit$rotation[, 1L], PC2 = fit$rotation[, 2L],
                         row.names = NULL)
  explained <- data.frame(PC = seq_along(fraction), fraction = fraction)
  list(scores = scores, loadings = loadings, explained = explained,
       dropped = features[!keep])
}

weighted_correlation <- function(X, weights, method = "pearson") {
  X <- as.matrix(X)
  if (!is.numeric(X) || any(!is.finite(X)) ||
      length(weights) != nrow(X) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("Invalid data or weights for correlation.", call. = FALSE)
  }
  if (!method %in% c("pearson", "spearman")) {
    stop("Unsupported correlation method.", call. = FALSE)
  }

  if (identical(method, "spearman")) {
    X <- apply(X, 2L, rank, ties.method = "average")
    X <- matrix(X, nrow = length(weights), dimnames = dimnames(X))
  }
  w <- as.numeric(weights) / sum(weights)
  means <- colSums(sweep(X, 1L, w, `*`))
  centered <- sweep(X, 2L, means, `-`)
  covariance <- crossprod(sweep(centered, 1L, sqrt(w), `*`))
  sds <- sqrt(pmax(diag(covariance), 0))
  denom <- outer(sds, sds)
  result <- covariance / denom
  result[denom < 1e-12] <- NA_real_
  # Assign back into the existing matrix: pmin/pmax with scalar first
  # arguments may discard the matrix dimensions and break dimnames().
  result[] <- pmin(1, pmax(-1, result))
  dimnames(result) <- list(colnames(X), colnames(X))
  result
}

run_correlation <- function(data, method = "pearson", weighted = FALSE,
                            focus = 2500, bandwidth = 750) {
  X <- as.matrix(data[, feature_names(data), drop = FALSE])
  if (!method %in% c("pearson", "spearman")) {
    stop("Unsupported correlation method.", call. = FALSE)
  }
  if (isTRUE(weighted)) {
    if (length(focus) != 1L || length(bandwidth) != 1L ||
        !is.finite(focus) || !is.finite(bandwidth) || bandwidth <= 0) {
      stop("Invalid Gaussian weighting settings.", call. = FALSE)
    }
    weights <- exp(-0.5 * ((data$time - focus) / bandwidth)^2)
    # Keep strictly positive weights even for far-away windows.
    weights <- pmax(weights, 1e-12)
  } else {
    weights <- rep(1, nrow(X))
  }
  mat <- weighted_correlation(X, weights, method)
  upper <- mat[upper.tri(mat)]
  finite_upper <- upper[is.finite(upper)]
  list(matrix = mat, effective_n = sum(weights)^2 / sum(weights^2),
       mean_abs = if (length(finite_upper)) mean(abs(finite_upper)) else NA_real_,
       method = method, weighted = isTRUE(weighted))
}

safe_result <- function(expr) {
  tryCatch(expr, error = function(e) {
    structure(list(message = conditionMessage(e)), class = "lab_error")
  })
}
