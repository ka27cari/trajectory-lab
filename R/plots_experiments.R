# Extended plots: isolated from server logic to keep the service testable.

plot_progressive_pca <- function(result) {
  current <- result$current
  reference <- result$reference
  ggplot2::ggplot() +
    ggplot2::geom_path(data = reference, ggplot2::aes(PC1, PC2),
                       colour = "#BBC8DA", linewidth = 1, linetype = "22") +
    ggplot2::geom_point(data = reference, ggplot2::aes(PC1, PC2),
                        colour = "#CBD5E3", size = 1.5, alpha = 0.5) +
    ggplot2::geom_path(data = current, ggplot2::aes(PC1, PC2),
                       colour = LAB_NAVY, linewidth = 0.95, alpha = 0.7) +
    ggplot2::geom_point(data = current, ggplot2::aes(PC1, PC2, colour = time),
                        size = 2.5, alpha = 0.92) +
    lab_time_scale() +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = sprintf("%d of %d features included", result$step, length(result$sequence)),
      subtitle = if (identical(result$mode, "refit"))
        "Dark = refitted PCA (aligned); dashed = all-feature reference" else
        "Dark = partial contribution in fixed axes; dashed = full scores",
      x = "Reference PC1 direction", y = "Reference PC2 direction") +
    theme_lab() + ggplot2::theme(legend.position = "bottom")
}

plot_loading_walk <- function(walk) {
  last <- tail(walk, 1L)
  d <- walk[-1L, , drop = FALSE]
  # Label a small readable subset rather than piling 40 labels on each other.
  every <- max(1L, ceiling(nrow(d) / 8L))
  labels <- d[d$step %% every == 0L | d$step == max(d$step), , drop = FALSE]
  ggplot2::ggplot(walk, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "#E5EBF4") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "#E5EBF4") +
    ggplot2::geom_path(colour = LAB_BLUE, linewidth = 0.85,
                       arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")) +
    ggplot2::geom_point(size = 2.1, colour = LAB_NAVY) +
    ggplot2::geom_point(data = last, size = 4, colour = LAB_CORAL) +
    ggplot2::geom_text(data = labels, ggplot2::aes(label = step),
                       size = 3, colour = LAB_NAVY, nudge_y = 0.04,
                       check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Cumulative loading walk",
                  subtitle = "Each arrow adds one feature loading in the FULL PCA basis",
                  x = "Cumulative PC1 loading", y = "Cumulative PC2 loading") +
    theme_lab()
}

plot_window_change_map <- function(summary, selected) {
  valid <- summary[is.finite(summary$variance_ratio) &
                     summary$variance_ratio >= 0, , drop = FALSE]
  if (!nrow(valid)) stop("No finite variance ratios are available.", call. = FALSE)
  # Cap log2 extreme ratios for legible axes, including exact zero late variance.
  valid$log_ratio <- pmax(-6, pmin(6, log2(pmax(valid$variance_ratio, 2^-6))))
  valid$highlight <- ifelse(valid$feature == selected, "Selected", "Other")
  ggplot2::ggplot(valid, ggplot2::aes(delta_mean, log_ratio)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B8C8D8", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, colour = "#B8C8D8", linetype = "dashed") +
    ggplot2::geom_point(data = valid[valid$highlight == "Other", , drop = FALSE],
                        colour = "#98A9BF", size = 2.4, alpha = 0.8) +
    ggplot2::geom_point(data = valid[valid$highlight == "Selected", , drop = FALSE],
                        colour = LAB_CORAL, size = 3.8) +
    ggplot2::labs(title = "Early vs late: mean and variance",
                  subtitle = "Below zero = lower late variance; ratios beyond 2^±6 are capped",
                  x = "Change in mean (late − early)",
                  y = "log2 variance ratio (late / early)") +
    theme_lab()
}

plot_half_distributions <- function(data, feature) {
  half <- floor(nrow(data) / 2L)
  d <- data.frame(
    period = factor(c(rep("Early", half), rep("Late", nrow(data) - half)),
                    levels = c("Early", "Late")),
    signal = data[[feature]])
  ggplot2::ggplot(d, ggplot2::aes(period, signal, fill = period)) +
    ggplot2::geom_boxplot(width = 0.42, alpha = 0.22, outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0,
                                                           seed = 11),
                        alpha = 0.32, size = 1.2) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 3,
                          fill = "white", colour = LAB_NAVY) +
    ggplot2::scale_fill_manual(values = c(Early = LAB_BLUE, Late = LAB_CORAL),
                               guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = feature, subtitle = "Diamonds mark the period means",
                  x = NULL, y = "Beta-like signal") + theme_lab()
}

plot_perturbed_trajectory <- function(paired, feature) {
  d <- data.frame(time = rep(paired$time, 2L),
                  signal = c(paired$baseline, paired$altered),
                  series = factor(rep(c("Original", "Experiment"),
                                      each = nrow(paired)),
                                  levels = c("Original", "Experiment")))
  ggplot2::ggplot(d, ggplot2::aes(time, signal, colour = series)) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::scale_colour_manual(values = c(Original = LAB_MUTED,
                                            Experiment = LAB_CORAL), name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = feature, subtitle = "Original and experimentally transformed signal",
                  x = "Time", y = "Beta-like signal") +
    theme_lab() + ggplot2::theme(legend.position = "bottom")
}

plot_perturbed_spread <- function(paired) {
  d <- data.frame(series = factor(rep(c("Original", "Experiment"),
                                      each = nrow(paired)),
                                  levels = c("Original", "Experiment")),
                  signal = c(paired$baseline, paired$altered))
  ggplot2::ggplot(d, ggplot2::aes(series, signal, fill = series)) +
    ggplot2::geom_boxplot(width = 0.5, alpha = 0.24, outlier.shape = NA) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.09, height = 0,
                                                           seed = 12),
                        alpha = 0.29, size = 1.2, colour = LAB_NAVY) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23,
                          size = 3.3, fill = "white", colour = LAB_NAVY) +
    ggplot2::scale_fill_manual(values = c(Original = LAB_MUTED,
                                          Experiment = LAB_CORAL), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Distribution before / after",
                  subtitle = "Mean marker = diamond; box width = IQR",
                  x = NULL, y = "Beta-like signal") + theme_lab()
}
