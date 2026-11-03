# Distance-based 2D ordinations and comparable descriptive diagnostics.
# Pure R calculations (MASS is a recommended R package, used for NMDS/Sammon).
# Rows of all working matrices are mapped objects; columns are measurements.

ordination_methods <- c("PCA", "PCoA", "NMDS", "Sammon")

distance_labels <- c(euclidean = "Euclidean", manhattan = "Manhattan",
                     canberra = "Canberra", bray = "Bray-Curtis",
                     correlation = "1 - Pearson correlation")

prepare_ordination <- function(data, orientation = "time", standardize = TRUE,
                               max_points = 120L) {
  if (!is.data.frame(data) || !"time" %in% names(data) ||
      !is.numeric(data$time) || any(!is.finite(data$time))) {
    stop("A numeric, finite time column is required.", call. = FALSE)
  }
  if (!orientation %in% c("time", "feature") ||
      length(max_points) != 1L || !is.finite(max_points) ||
      max_points < 10L || max_points > 200L ||
      max_points != as.integer(max_points) ||
      length(standardize) != 1L || is.na(standardize)) {
    stop("Invalid ordination orientation or point limit.", call. = FALSE)
  }
  vars <- feature_names(data)
  if (length(vars) < 2L) stop("At least two feature columns are required.", call. = FALSE)
  X <- as.matrix(data[, vars, drop = FALSE])
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("Ordination requires finite numeric measurements.", call. = FALSE)
  }
  # A constant feature adds no pairwise distance and cannot be standardized.
  keep <- apply(X, 2L, stats::sd) > 1e-10
  removed <- colnames(X)[!keep]
  X <- X[, keep, drop = FALSE]
  if (ncol(X) < 2L) {
    stop("At least two non-constant features are needed.", call. = FALSE)
  }
  if (isTRUE(standardize)) {
    X <- scale(X, center = TRUE, scale = TRUE)
  }
  if (identical(orientation, "time")) {
    idx <- if (nrow(X) > max_points)
      unique(as.integer(round(seq(1, nrow(X), length.out = max_points)))) else
      seq_len(nrow(X))
    Y <- X[idx, , drop = FALSE]
    labels <- sprintf("T%03d", idx)
    time <- data$time[idx]
    selected <- length(idx)
  } else {
    Y <- t(X)
    labels <- rownames(Y)
    time <- rep(NA_real_, nrow(Y))
    selected <- nrow(Y)
  }
  if (nrow(Y) < 3L || ncol(Y) < 2L) {
    stop("Two-dimensional maps require at least 3 objects and 2 measurements.", call. = FALSE)
  }
  rownames(Y) <- labels
  list(X = unname(Y), labels = labels, time = time,
       orientation = orientation, scaled = isTRUE(standardize),
       n_objects = nrow(Y), total_objects = if (orientation == "time") nrow(data) else length(vars),
       features_removed = removed,
       n_measurements = ncol(Y), max_points = max_points)
}

make_dissimilarity <- function(prepared, method = "euclidean") {
  if (length(method) != 1L || !method %in% names(distance_labels)) {
    stop("Unknown distance metric.", call. = FALSE)
  }
  X <- prepared$X
  if (identical(method, "bray")) {
    if (any(X < 0)) {
      stop("Bray-Curtis requires nonnegative values: switch Standardize to No.", call. = FALSE)
    }
    n <- nrow(X)
    d <- matrix(0, n, n)
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        denom <- sum(X[i, ] + X[j, ])
        if (denom <= 0) {
          stop("Bray-Curtis is undefined for two all-zero objects.", call. = FALSE)
        }
        d[i, j] <- d[j, i] <- sum(abs(X[i, ] - X[j, ])) / denom
      }
    }
    D <- stats::as.dist(d)
  } else if (identical(method, "correlation")) {
    sd_rows <- apply(X, 1L, stats::sd)
    if (any(!is.finite(sd_rows)) || any(sd_rows < 1e-12)) {
      stop("Correlation distance is undefined for an object with constant measurements.",
           call. = FALSE)
    }
    d <- 1 - stats::cor(t(X), method = "pearson")
    d[] <- pmax(0, pmin(2, d))
    diag(d) <- 0
    D <- stats::as.dist(d)
  } else {
    D <- stats::dist(X, method = method)
  }
  values <- as.numeric(D)
  if (!length(values) || any(!is.finite(values)) || any(values < 0) ||
      !any(values > 1e-12)) {
    stop("This selection has no usable pairwise differences.", call. = FALSE)
  }
  attr(D, "Labels") <- prepared$labels
  D
}

# Distances have different scales. 2D fit errors are computed after an optimal
# nonnegative scalar adjustment, so only *relative distance distortion* remains.
ordination_quality <- function(original, points) {
  target <- as.numeric(original)
  mapped <- as.numeric(stats::dist(points))
  denom <- sum(mapped^2)
  factor <- if (denom <= 1e-15) 0 else sum(target * mapped) / denom
  relative_error <- sqrt(sum((target - factor * mapped)^2) / sum(target^2))
  rank <- if (stats::sd(mapped) <= 1e-14 || stats::sd(target) <= 1e-14)
    NA_real_ else suppressWarnings(stats::cor(target, mapped, method = "spearman"))
  list(distance_error = relative_error, distance_rank = rank,
       scale_factor = factor)
}

# One method may be impossible (e.g. tied distances for NMDS) without blocking
# the other maps. Return the individual error as ordinary status text.
fit_ordination <- function(prepared, dissimilarity, method) {
  if (!method %in% ordination_methods) stop("Unknown ordination method.", call. = FALSE)
  X <- prepared$X
  n <- nrow(X)
  extra <- ""
  if (identical(method, "PCA")) {
    fit <- stats::prcomp(X, center = TRUE, scale. = FALSE)
    points <- fit$x[, 1:2, drop = FALSE]
    explain <- (fit$sdev[1:2]^2) / sum(fit$sdev^2)
    extra <- sprintf("PC1 %.1f%% / PC2 %.1f%%", 100 * explain[1L], 100 * explain[2L])
  } else if (identical(method, "PCoA")) {
    fit <- suppressWarnings(stats::cmdscale(dissimilarity, k = 2L, eig = TRUE))
    points <- as.matrix(fit$points)
    if (is.null(dim(points)) || ncol(points) < 2L) {
      stop("PCoA has fewer than two positive axes for these distances.", call. = FALSE)
    }
    eig <- fit$eig
    positive <- sum(eig[eig > 0])
    negative <- abs(sum(eig[eig < 0]))
    share <- if (positive > 0) sum(head(eig[eig > 0], 2L)) / positive else 0
    extra <- sprintf("Positive inertia %.1f%%; negative eigenvalue mass %.1f%%",
                     100 * share, 100 * negative / max(positive + negative, 1e-12))
  } else {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Install the recommended R package MASS for NMDS and Sammon.", call. = FALSE)
    }
    if (any(as.numeric(dissimilarity) <= 1e-12)) {
      stop(paste(method, "needs strictly positive off-diagonal distances; duplicate objects were found."),
           call. = FALSE)
    }
    initial <- suppressWarnings(stats::cmdscale(dissimilarity, k = 2L))
    if (!is.matrix(initial) || ncol(initial) != 2L ||
        anyDuplicated(as.data.frame(signif(initial, 12)))) {
      stop(paste(method, "needs a non-degenerate 2D initial map."), call. = FALSE)
    }
    if (identical(method, "NMDS")) {
      fit <- suppressWarnings(MASS::isoMDS(dissimilarity, y = initial, k = 2L,
                                           maxit = 80L, trace = FALSE))
      points <- fit$points
      extra <- sprintf("Kruskal stress %.2f%% (method-specific)", fit$stress)
    } else {
      fit <- suppressWarnings(MASS::sammon(dissimilarity, y = initial, k = 2L,
                                           niter = 100L, trace = FALSE))
      points <- fit$points
      extra <- sprintf("Sammon stress %.4f (method-specific)", fit$stress)
    }
  }
  if (!is.matrix(points) || any(dim(points) != c(n, 2L)) || any(!is.finite(points))) {
    stop("The method did not produce valid two-dimensional coordinates.", call. = FALSE)
  }
  colnames(points) <- c("Dim1", "Dim2")
  qual <- ordination_quality(dissimilarity, points)
  mapped <- data.frame(label = prepared$labels,
                       time = prepared$time,
                       Dim1 = points[, 1L], Dim2 = points[, 2L],
                       method = method, stringsAsFactors = FALSE)
  list(scores = mapped, quality = qual, detail = extra)
}

run_ordination_suite <- function(data, distance = "euclidean", standardize = TRUE,
                                 orientation = "time", max_points = 120L) {
  prepared <- prepare_ordination(data, orientation, standardize, max_points)
  diss <- make_dissimilarity(prepared, distance)
  fits <- stats::setNames(lapply(ordination_methods, function(method) {
    tryCatch(fit_ordination(prepared, diss, method),
             error = function(e) list(error = conditionMessage(e)))
  }), ordination_methods)
  results <- lapply(ordination_methods, function(method) {
    obj <- fits[[method]]
    success <- is.null(obj$error)
    data.frame(method = method,
               distance_rank = if (success) obj$quality$distance_rank else NA_real_,
               distance_error = if (success) obj$quality$distance_error else NA_real_,
               detail = if (success) obj$detail else obj$error,
               status = if (success) "OK" else "Unavailable",
               stringsAsFactors = FALSE)
  })
  list(data = prepared, distance = diss, distance_method = distance,
       fits = fits, summary = do.call(rbind, results))
}

ordination_scores <- function(suite) {
  valid <- suite$fits[vapply(suite$fits, function(x) is.null(x$error), logical(1))]
  if (!length(valid)) stop("No ordination methods produced a usable map.", call. = FALSE)
  do.call(rbind, lapply(valid, function(x) x$scores))
}

ordination_pair_data <- function(suite, method, n_max = 2500L) {
  fit <- suite$fits[[method]]
  if (is.null(fit) || !is.null(fit$error)) stop("No coordinates for that method.", call. = FALSE)
  original <- as.numeric(suite$distance)
  mapped <- as.numeric(stats::dist(fit$scores[, c("Dim1", "Dim2")]))
  idx <- if (length(original) > n_max) unique(as.integer(round(
    seq(1, length(original), length.out = n_max)))) else seq_along(original)
  data.frame(dissimilarity = original[idx],
             embedded_distance = mapped[idx] * fit$quality$scale_factor)
}

ordination_distance_preview <- function(suite, limit = 48L) {
  D <- as.matrix(suite$distance)
  idx <- if (nrow(D) > limit) unique(as.integer(round(seq(1, nrow(D),
                                                       length.out = limit)))) else seq_len(nrow(D))
  D[idx, idx, drop = FALSE]
}
