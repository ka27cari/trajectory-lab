# Synthetic trajectories are illustrative beta-like signals, not empirical CpG data.

with_local_seed <- function(seed, expression) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed) || seed < 1L) {
    stop("The random seed must be a positive integer.", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expression)
}

simulate_one <- function(time, midpoint, signed_amplitude, baseline,
                         transition_days, noise_sd, seed) {
  stopifnot(length(time) >= 2L, all(is.finite(time)))
  if (!is.finite(transition_days) || transition_days <= 0 ||
      !is.finite(noise_sd) || noise_sd < 0) {
    stop("Transition width must be positive, and noise cannot be negative.", call. = FALSE)
  }

  with_local_seed(seed, {
    fraction <- stats::plogis(4 * (time - midpoint) / transition_days)
    signal <- baseline + signed_amplitude * (fraction - 0.5)
    # Clipping is part of this illustrative generator and can affect correlation.
    pmin(1, pmax(0, signal + stats::rnorm(length(time), sd = noise_sd)))
  })
}

simulate_dataset <- function(n_features = 12L, n_timepoints = 140L,
                             pattern = "mixed", amplitude = 0.52,
                             transition_days = 900, noise_sd = 0.035,
                             seed = 42L) {
  n_features <- as.integer(n_features)
  n_timepoints <- as.integer(n_timepoints)
  seed <- as.integer(seed)
  if (anyNA(c(n_features, n_timepoints, seed, amplitude,
              transition_days, noise_sd)) ||
      !is.finite(amplitude) || amplitude <= 0 || amplitude > 0.9 ||
      !is.finite(transition_days) || transition_days <= 0 ||
      !is.finite(noise_sd) || noise_sd < 0 || noise_sd > 0.25 ||
      n_features < 2L || n_features > 30L ||
      n_timepoints < 20L || n_timepoints > 500L ||
      seed < 1L || seed > 1000000L ||
      !pattern %in% c("mixed", "early", "late", "synchronous")) {
    stop("Invalid simulation settings.", call. = FALSE)
  }

  config <- list(n_features = n_features, n_timepoints = n_timepoints,
                 pattern = pattern, amplitude = amplitude,
                 transition_days = transition_days, noise_sd = noise_sd,
                 seed = seed)
  time <- seq(0, 5000, length.out = n_timepoints)

  parameters <- with_local_seed(seed, {
    centers <- switch(pattern,
      mixed = stats::runif(n_features, 500, 4500),
      early = stats::runif(n_features, 600, 1950),
      late = stats::runif(n_features, 3050, 4400),
      synchronous = stats::runif(n_features, 2180, 2820)
    )
    directions <- sample(c(-1, 1), size = n_features, replace = TRUE)
    data.frame(
      feature = sprintf("CpG_%02d", seq_len(n_features)),
      midpoint = round(centers, 1),
      signed_amplitude = round(directions * amplitude * stats::runif(n_features, 0.65, 1), 4),
      baseline = round(stats::runif(n_features, 0.43, 0.57), 4),
      transition_days = round(transition_days * stats::runif(n_features, 0.8, 1.2), 1),
      noise_sd = rep(noise_sd, n_features),
      noise_seed = seed + seq_len(n_features),
      source = rep("simulated", n_features),
      stringsAsFactors = FALSE
    )
  })

  curves <- vapply(seq_len(n_features), function(i) {
    simulate_one(
      time, parameters$midpoint[i], parameters$signed_amplitude[i],
      parameters$baseline[i], parameters$transition_days[i],
      parameters$noise_sd[i], parameters$noise_seed[i]
    )
  }, numeric(length(time)))

  colnames(curves) <- parameters$feature
  data <- data.frame(time = time, curves, check.names = FALSE)
  list(data = data, metadata = parameters, config = config, source = "simulated")
}

add_custom_feature <- function(dataset, midpoint, signed_amplitude,
                               baseline, transition_days, noise_sd) {
  stopifnot(identical(dataset$source, "simulated"))
  if (ncol(dataset$data) >= 41L) {
    stop("A maximum of 40 features is supported.", call. = FALSE)
  }
  params <- c(midpoint, signed_amplitude, baseline, transition_days, noise_sd)
  if (any(!is.finite(params)) || midpoint < 0 || midpoint > 5000 ||
      abs(signed_amplitude) > 0.9 || baseline < 0 || baseline > 1 ||
      transition_days < 50 || noise_sd < 0 || noise_sd > 0.25) {
    stop("Invalid custom feature settings.", call. = FALSE)
  }

  feature <- sprintf("CpG_%02d", ncol(dataset$data))
  noise_seed <- dataset$config$seed + 10000L + ncol(dataset$data)
  dataset$data[[feature]] <- simulate_one(
    dataset$data$time, midpoint, signed_amplitude, baseline,
    transition_days, noise_sd, noise_seed
  )
  dataset$metadata <- rbind(dataset$metadata, data.frame(
    feature = feature, midpoint = midpoint,
    signed_amplitude = signed_amplitude, baseline = baseline,
    transition_days = transition_days, noise_sd = noise_sd,
    noise_seed = noise_seed, source = "custom", stringsAsFactors = FALSE
  ))
  dataset
}

undo_custom_feature <- function(dataset) {
  stopifnot(identical(dataset$source, "simulated"))
  if (nrow(dataset$metadata) <= dataset$config$n_features) {
    stop("No custom feature to remove.", call. = FALSE)
  }
  last_feature <- tail(dataset$metadata$feature, 1L)
  dataset$data[[last_feature]] <- NULL
  dataset$metadata <- head(dataset$metadata, -1L)
  dataset
}
