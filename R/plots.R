LAB_BLUE <- "#4777E8"
LAB_CORAL <- "#EA876B"
LAB_NAVY <- "#152941"
LAB_MUTED <- "#8293AD"

theme_lab <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = LAB_NAVY),
      plot.title = ggplot2::element_text(face = "bold", size = 14, margin = ggplot2::margin(b = 7)),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#61738E"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#EBEFF6", linewidth = 0.35),
      axis.title = ggplot2::element_text(size = 10, colour = "#50627D"),
      axis.text = ggplot2::element_text(colour = "#63748D", size = 9),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 9),
      legend.key.height = grid::unit(0.32, "cm"),
      legend.spacing.x = grid::unit(0.16, "cm"),
      plot.margin = ggplot2::margin(14, 18, 14, 12)
    )
}

# Keep color-bar ticks and labels legible in narrow cards.
lab_time_scale <- function() {
  ggplot2::scale_colour_gradient(
    low = LAB_BLUE, high = LAB_CORAL, name = "Time",
    breaks = function(limits) pretty(limits, n = 3),
    labels = function(x) ifelse(abs(x) >= 1000,
                                  paste0(format(round(x / 1000, 1), trim = TRUE), "k"),
                                  format(signif(x, 3), trim = TRUE)),
    guide = ggplot2::guide_colourbar(
      title.position = "top", title.hjust = 0,
      barwidth = grid::unit(4.1, "cm"),
      barheight = grid::unit(0.24, "cm")))
}

long_profiles <- function(data, features = feature_names(data)) {
  n <- nrow(data)
  data.frame(
    time = rep(data$time, times = length(features)),
    feature = rep(features, each = n),
    beta = unlist(data[, features, drop = FALSE], use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

plot_trajectories <- function(data, selected) {
  long <- long_profiles(data)
  long$emphasis <- ifelse(long$feature == selected, "Selected", "Other features")
  ggplot2::ggplot(long, ggplot2::aes(x = time, y = beta, group = feature)) +
    ggplot2::geom_line(data = long[long$emphasis == "Other features", , drop = FALSE],
                       colour = LAB_MUTED, linewidth = 0.65, alpha = 0.33) +
    ggplot2::geom_line(data = long[long$emphasis == "Selected", , drop = FALSE],
                       colour = LAB_BLUE, linewidth = 1.4) +
    ggplot2::labs(title = "Trajectory overview",
                  subtitle = paste("Highlighted feature:", selected),
                  x = "Time (arbitrary days)", y = "Beta-like signal") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_lab()
}

plot_feature <- function(data, feature, metadata_row = NULL) {
  selected <- data.frame(time = data$time, beta = data[[feature]])
  p <- ggplot2::ggplot(selected, ggplot2::aes(time, beta)) +
    ggplot2::geom_area(fill = LAB_BLUE, alpha = 0.09) +
    ggplot2::geom_line(colour = LAB_BLUE, linewidth = 1.05) +
    ggplot2::geom_point(colour = LAB_BLUE, size = 0.7, alpha = 0.55) +
    ggplot2::labs(title = feature, subtitle = "Observed points in the selected window",
                  x = "Time", y = "Beta-like signal") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) + theme_lab()
  if (!is.null(metadata_row) && nrow(metadata_row) > 0L &&
      is.finite(metadata_row$midpoint[1L])) {
    p <- p + ggplot2::geom_vline(xintercept = metadata_row$midpoint[1L],
                                 linetype = "dashed", colour = LAB_CORAL, linewidth = 0.8)
  }
  p
}

plot_pca_scores <- function(result) {
  d <- result$scores
  explained <- result$explained$fraction
  ggplot2::ggplot(d, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_path(colour = "#93A5BC", linewidth = 0.85, alpha = 0.65) +
    ggplot2::geom_point(ggplot2::aes(colour = time), size = 2.6, alpha = 0.92) +
    lab_time_scale() +
    ggplot2::labs(
      title = "PCA score trajectory",
      subtitle = "Each point is a time point; each feature is a variable",
      x = sprintf("PC1 (%.1f%%)", 100 * explained[1L]),
      y = sprintf("PC2 (%.1f%%)", 100 * explained[2L])
    ) + theme_lab() +
    ggplot2::theme(legend.position = "bottom")
}

plot_pca_scree <- function(result) {
  d <- head(result$explained, 12L)
  ggplot2::ggplot(d, ggplot2::aes(x = factor(PC), y = fraction)) +
    ggplot2::geom_col(fill = LAB_BLUE, width = 0.62, alpha = 0.84) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                                limits = c(0, max(d$fraction, na.rm = TRUE) * 1.16)) +
    ggplot2::scale_x_discrete(breaks = as.character(seq.int(1L, nrow(d),
                                           by = if (nrow(d) > 8L) 2L else 1L))) +
    ggplot2::labs(title = "Explained variance", subtitle = "Per principal component",
                  x = "Component", y = "Variance share") +
    theme_lab()
}

plot_pca_loadings <- function(result) {
  d <- result$loadings
  d <- d[order(abs(d$PC1), decreasing = TRUE), , drop = FALSE]
  d <- head(d, 12L)
  d$feature <- stats::reorder(d$feature, d$PC1)
  ggplot2::ggplot(d, ggplot2::aes(feature, PC1, fill = PC1 > 0)) +
    ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("TRUE" = LAB_CORAL, "FALSE" = LAB_BLUE)) +
    ggplot2::labs(title = "PC1 feature loadings",
                  subtitle = "Signed coefficients, not time-point coordinates",
                  x = NULL, y = "Loading") + theme_lab()
}

plot_correlation_heatmap <- function(result) {
  mat <- result$matrix
  features <- rownames(mat)
  d <- expand.grid(x = features, y = features, stringsAsFactors = FALSE)
  d$r <- as.vector(mat)
  d$x <- factor(d$x, levels = features)
  d$y <- factor(d$y, levels = rev(features))

  ggplot2::ggplot(d, ggplot2::aes(x, y, fill = r)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(low = LAB_BLUE, mid = "#FFFFFF",
                                  high = LAB_CORAL, midpoint = 0,
                                  limits = c(-1, 1), breaks = c(-1, 0, 1),
                                  na.value = "#E8EAF0", name = "Correlation",
                                  guide = ggplot2::guide_colourbar(
                                    title.position = "top", title.hjust = 0,
                                    barwidth = grid::unit(4.1, "cm"),
                                    barheight = grid::unit(0.24, "cm"))) +
    ggplot2::scale_x_discrete(
      breaks = features[seq.int(1L, length(features),
                                by = max(1L, ceiling(length(features) / 8L)))],
      drop = FALSE) +
    ggplot2::scale_y_discrete(
      breaks = rev(features)[seq.int(1L, length(features),
                                     by = max(1L, ceiling(length(features) / 12L)))],
      drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Feature correlation map",
                  subtitle = "Click a tile to inspect a pair",
                  x = NULL, y = NULL) +
    theme_lab() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   axis.text.y = ggplot2::element_text(size = 8),
                   panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom")
}

plot_pair_scatter <- function(data, feature_a, feature_b) {
  d <- data.frame(time = data$time, A = data[[feature_a]], B = data[[feature_b]])
  ggplot2::ggplot(d, ggplot2::aes(A, B)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         colour = "#CCD5E3", linewidth = 0.55, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = time), size = 2, alpha = 0.8) +
    lab_time_scale() +
    ggplot2::labs(title = "Pairwise scatter", x = feature_a, y = feature_b) +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_lab() + ggplot2::theme(legend.position = "bottom")
}

plot_pair_trajectories <- function(data, feature_a, feature_b) {
  d <- long_profiles(data, features = unique(c(feature_a, feature_b)))
  ggplot2::ggplot(d, ggplot2::aes(time, beta, colour = feature)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c(LAB_BLUE, LAB_CORAL)[seq_along(unique(d$feature))], unique(d$feature))) +
    ggplot2::labs(title = "Selected pair over time", x = "Time", y = "Beta-like signal",
                  colour = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) + theme_lab() +
    ggplot2::theme(legend.position = "bottom")
}
