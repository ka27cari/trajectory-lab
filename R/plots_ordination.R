# Visualizations for one-screen ordination comparison. Uses same theme_lab
# and time color scale as all other charts.

plot_ordination_comparison <- function(suite) {
  d <- ordination_scores(suite)
  order <- suite$summary$method[suite$summary$status == "OK"]
  d$method <- factor(d$method, levels = order)
  subtitle <- sprintf("%s distance  ·  %s  ·  %d/%d objects shown  ·  axes are method-specific",
                      distance_labels[[suite$distance_method]],
                      if (suite$data$scaled) "standardized" else "raw values",
                      suite$data$n_objects, suite$data$total_objects)
  if (identical(suite$data$orientation, "time")) {
    p <- ggplot2::ggplot(d, ggplot2::aes(Dim1, Dim2)) +
      ggplot2::geom_path(ggplot2::aes(group = method), colour = "#9BAFC4",
                         linewidth = 0.55, alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(colour = time), size = 1.85, alpha = 0.9) +
      lab_time_scale()
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(Dim1, Dim2)) +
      ggplot2::geom_point(colour = LAB_BLUE, size = 2.6, alpha = 0.9) +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 2.65,
                         nudge_y = 0.06, check_overlap = TRUE, colour = LAB_NAVY)
  }
  p + ggplot2::facet_wrap(~method, ncol = 2L, scales = "free") +
    ggplot2::labs(title = "Four ways to map the same objects", subtitle = subtitle,
                  x = "Dimension 1", y = "Dimension 2") +
    ggplot2::scale_x_continuous(breaks = function(x) pretty(x, n = 3)) +
    ggplot2::scale_y_continuous(breaks = function(x) pretty(x, n = 3)) +
    theme_lab() +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold", size = 11,
                                                       margin = ggplot2::margin(t = 8, b = 8)),
                   strip.background = ggplot2::element_rect(fill = "#EDF3FF", colour = NA),
                   panel.spacing = grid::unit(0.9, "lines"),
                   axis.text = ggplot2::element_text(size = 8),
                   plot.subtitle = ggplot2::element_text(size = 9))
}

plot_ordination_focus <- function(suite, method) {
  fit <- suite$fits[[method]]
  if (is.null(fit) || !is.null(fit$error)) stop("Select a successful method.", call. = FALSE)
  d <- fit$scores
  if (identical(suite$data$orientation, "time")) {
    p <- ggplot2::ggplot(d, ggplot2::aes(Dim1, Dim2)) +
      ggplot2::geom_path(colour = LAB_MUTED, linewidth = 0.85, alpha = 0.7) +
      ggplot2::geom_point(ggplot2::aes(colour = time), size = 2.6, alpha = 0.92) +
      lab_time_scale()
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(Dim1, Dim2)) +
      ggplot2::geom_point(size = 3, colour = LAB_BLUE) +
      ggplot2::geom_text(ggplot2::aes(label = label),
                         nudge_y = 0.04, check_overlap = TRUE, size = 3)
  }
  p + ggplot2::labs(title = paste(method, "2D map"), subtitle = fit$detail,
                    x = "Dimension 1", y = "Dimension 2") +
    theme_lab() + ggplot2::theme(legend.position = "bottom")
}

plot_ordination_shepard <- function(suite, method) {
  d <- ordination_pair_data(suite, method)
  fit <- suite$fits[[method]]
  ggplot2::ggplot(d, ggplot2::aes(dissimilarity, embedded_distance)) +
    ggplot2::geom_point(colour = LAB_BLUE, size = 1, alpha = 0.24) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = LAB_CORAL,
                         linetype = "dashed", linewidth = 0.8) +
    ggplot2::labs(title = "Dissimilarity vs 2D distance",
      subtitle = sprintf("Rank %.3f  ·  relative distance error %.3f",
                        fit$quality$distance_rank, fit$quality$distance_error),
      x = "Requested pairwise dissimilarity",
      y = "2D Euclidean distance (scale adjusted)") +
    theme_lab()
}

plot_dissimilarity_matrix <- function(suite) {
  D <- ordination_distance_preview(suite)
  labels <- rownames(D)
  if (is.null(labels)) labels <- as.character(seq_len(nrow(D)))
  d <- expand.grid(x = labels, y = labels, stringsAsFactors = FALSE)
  d$value <- as.vector(D)
  d$x <- factor(d$x, levels = labels)
  d$y <- factor(d$y, levels = rev(labels))
  ticks <- max(1L, ceiling(length(labels) / 8L))
  ggplot2::ggplot(d, ggplot2::aes(x, y, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_gradient(low = "#F1F6FF", high = LAB_BLUE,
                                 name = "Distance",
      guide = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0,
        barwidth = grid::unit(4.1, "cm"), barheight = grid::unit(0.23, "cm"))) +
    ggplot2::scale_x_discrete(breaks = labels[seq.int(1L, length(labels), by = ticks)], drop = FALSE) +
    ggplot2::scale_y_discrete(breaks = rev(labels)[seq.int(1L, length(labels), by = ticks)], drop = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Pairwise dissimilarities",
                  subtitle = sprintf("Showing %d of %d objects", nrow(D), suite$data$n_objects),
                  x = NULL, y = NULL) + theme_lab() +
    ggplot2::theme(legend.position = "bottom", panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                   axis.text.y = ggplot2::element_text(size = 7))
}
