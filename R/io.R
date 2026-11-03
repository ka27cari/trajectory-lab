read_beta_csv <- function(path) {
  if (!file.exists(path) || file.info(path)$size > 3 * 1024^2) {
    stop("Upload a CSV no larger than 3 MB.", call. = FALSE)
  }
  raw <- tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM"),
    error = function(e) stop("Unable to read CSV: ", conditionMessage(e), call. = FALSE)
  )
  if (!"time" %in% names(raw) || anyDuplicated(names(raw))) {
    stop("CSV requires a unique column named 'time'.", call. = FALSE)
  }
  features <- setdiff(names(raw), "time")
  if (length(features) < 2L || length(features) > 40L ||
      nrow(raw) < 5L || nrow(raw) > 1000L) {
    stop("CSV needs 5-1000 rows and 2-40 feature columns.", call. = FALSE)
  }
  if (!all(vapply(raw, is.numeric, logical(1L))) ||
      any(!is.finite(as.matrix(raw)))) {
    stop("All values must be finite numeric values without missing entries.", call. = FALSE)
  }
  if (any(diff(raw$time) <= 0)) {
    stop("'time' must be strictly increasing with no duplicates.", call. = FALSE)
  }
  X <- as.matrix(raw[, features, drop = FALSE])
  if (any(X < 0 | X > 1)) {
    stop("Feature values must be beta-scale signals between 0 and 1.", call. = FALSE)
  }
  metadata <- data.frame(
    feature = features, midpoint = NA_real_, signed_amplitude = NA_real_,
    baseline = NA_real_, transition_days = NA_real_, noise_sd = NA_real_,
    noise_seed = NA_integer_, source = "uploaded", stringsAsFactors = FALSE
  )
  list(data = raw[, c("time", features), drop = FALSE], metadata = metadata,
       config = NULL, source = "uploaded")
}

make_session_notes <- function(dataset, observation_window, method,
                               weighted, focus, bandwidth) {
  base <- c(
    "TRAJECTORY LAB | analysis settings",
    "=================================",
    "All profiles are illustrative unless you supplied your own CSV.",
    sprintf("Data source: %s", dataset$source),
    sprintf("Features: %d", length(feature_names(dataset$data))),
    sprintf("Time points (full dataset): %d", nrow(dataset$data)),
    sprintf("Selected window: %.6f to %.6f", observation_window[1], observation_window[2]),
    sprintf("Correlation: %s", method),
    sprintf("Gaussian time weighting: %s", if (weighted) "on" else "off")
  )
  if (weighted) {
    base <- c(base, sprintf("Weight center: %.6f", focus),
              sprintf("Weight bandwidth (SD): %.6f", bandwidth))
  }
  if (!is.null(dataset$config)) {
    base <- c(base, "", "Base simulator configuration:",
              capture.output(dput(dataset$config)))
  }
  c(base, "", "Download metadata.csv to capture each feature's generator parameters.",
    "Reproducibility also depends on the R version and RNG implementation.",
    "No biological associations or scientific validation are claimed.")
}
