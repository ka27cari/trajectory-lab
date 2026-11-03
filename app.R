# Trajectory Lab | a research-oriented synthetic time-course explorer.
# Launch from the repository root: shiny::runApp(".")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
})

source("R/simulation.R", local = TRUE)
source("R/analysis.R", local = TRUE)
source("R/experiments.R", local = TRUE)
source("R/ordination.R", local = TRUE)
source("R/io.R", local = TRUE)
source("R/plots.R", local = TRUE)
source("R/plots_experiments.R", local = TRUE)
source("R/plots_ordination.R", local = TRUE)

options(shiny.maxRequestSize = 3 * 1024^2)

metric_tile <- function(kicker, value, note) {
  tags$div(class = "metric",
           tags$div(class = "metric-kicker", kicker),
           tags$div(class = "metric-value", value),
           tags$div(class = "metric-note", note))
}

panel <- function(title, subtitle = NULL, ...) {
  tags$div(class = "lab-panel",
           tags$div(class = "panel-heading", tags$h3(title),
                    if (!is.null(subtitle)) tags$p(subtitle)),
           ...)
}

ui <- bslib::page_sidebar(
  title = "Trajectory Lab",
  theme = bslib::bs_theme(version = 5, bg = "#F4F7FC", fg = "#172941",
                          primary = "#4777E8"),
  fillable = FALSE,
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
    includeCSS("www/styles.css")
  ),
  sidebar = bslib::sidebar(
    width = 320,
    tags$div(class = "sidebar-intro",
             tags$div(class = "eyebrow", "CONTROL CENTER"),
             tags$p("Shape a synthetic signal, analyze its structure, and export what you see.")),
    radioButtons("source", "Data source", inline = TRUE,
                 choices = c("Simulator" = "simulated", "Upload CSV" = "uploaded"),
                 selected = "simulated"),
    conditionalPanel("input.source === 'simulated'",
      tags$div(class = "section-label", "01 / GENERATOR"),
      selectInput("pattern", "Transition pattern",
                  choices = c("Mixed timing" = "mixed", "Early transitions" = "early",
                              "Late transitions" = "late", "Synchronous transitions" = "synchronous"),
                  selected = "mixed"),
      sliderInput("n_features", "Number of features", min = 2, max = 30, value = 12, step = 1),
      sliderInput("n_timepoints", "Time points", min = 20, max = 300, value = 140, step = 10),
      sliderInput("amplitude", "Signal amplitude", min = 0.1, max = 0.9,
                  value = 0.52, step = 0.02),
      sliderInput("transition_days", "Transition width (days)", min = 150,
                  max = 2200, value = 900, step = 50),
      sliderInput("noise_sd", "Noise standard deviation", min = 0, max = 0.2,
                  value = 0.035, step = 0.005),
      numericInput("seed", "Random seed", value = 42, min = 1, max = 1000000, step = 1),
      actionButton("generate", "Regenerate dataset", class = "btn-primary full-button"),
      tags$details(class = "custom-feature",
        tags$summary("+ Add a custom feature"),
        sliderInput("custom_midpoint", "Midpoint", min = 0, max = 5000,
                    value = 2500, step = 100),
        sliderInput("custom_amplitude", "Signed amplitude", min = -0.9,
                    max = 0.9, value = -0.45, step = 0.05),
        sliderInput("custom_baseline", "Baseline", min = 0.1, max = 0.9,
                    value = 0.5, step = 0.05),
        sliderInput("custom_width", "Transition width", min = 100, max = 2200,
                    value = 850, step = 50),
        actionButton("add_feature", "Add feature", class = "btn-outline-primary"),
        actionButton("undo_feature", "Undo last", class = "btn-outline-secondary")
      )
    ),
    conditionalPanel("input.source === 'uploaded'",
      tags$div(class = "section-label", "01 / YOUR DATA"),
      fileInput("user_csv", "Upload beta-scale CSV", accept = c(".csv", "text/csv")),
      tags$p(class = "hint", "Required: a numeric 'time' column, followed by 2–40 feature columns (0–1 values). Five to 1,000 ordered time points; 3 MB max."),
      tags$p(class = "hint", "Uploaded files are analyzed in the current Shiny session and are not saved by the app code.")
    ),
    tags$div(class = "section-label", "02 / ANALYSIS"),
    sliderInput("time_window", "Observation window", min = 0, max = 5000,
                value = c(0, 5000), step = 50),
    selectInput("cor_method", "Correlation measure",
                choices = c("Pearson" = "pearson", "Spearman" = "spearman")),
    checkboxInput("weighted", "Gaussian time weighting", value = FALSE),
    conditionalPanel("input.weighted === true",
      numericInput("focus_time", "Weight center (time)", value = 2500, min = -1000000),
      numericInput("bandwidth", "Weight bandwidth (SD)", value = 750, min = 0.000001)
    ),
    tags$div(class = "sidebar-footnote", "Synthetic examples only. Not a diagnostic or validated inference tool.")
  ),
  tags$div(class = "hero",
    tags$div(class = "hero-copy",
      tags$div(class = "eyebrow", "INTERACTIVE RESEARCH SANDBOX · v1.2.0"),
      tags$h1("Explore the shape of change."),
      tags$p("Shape trajectories, follow PCA changes, test mean and spread, and compare four distance-aware maps in one view.")
    ),
    tags$div(class = "hero-mark", "TL", `aria-hidden` = "true")
  ),
  uiOutput("metrics"),
  bslib::navset_card_tab(
    id = "view",
    bslib::nav_panel("Explore",
      tags$div(class = "lab-grid two-one",
        panel("Signals", "Select any feature to inspect its time course.",
          selectInput("selected_feature", label = NULL, choices = NULL, width = "200px"),
          plotOutput("trajectory_plot", height = "410px")),
        panel("Feature inspector", "Parameter metadata and observed signal.",
          uiOutput("feature_stats"), plotOutput("feature_plot", height = "320px"))
      )
    ),
    bslib::nav_panel("PCA",
      tags$div(class = "lab-grid two-one",
        panel("Principal components", "Scaled feature values in the selected time window.",
          plotOutput("pca_plot", height = "440px")),
        panel("Variance decomposition", "Percentage of variance described by each axis.",
          plotOutput("scree_plot", height = "440px"))
      ),
      tags$div(class = "lab-grid one",
        panel("Feature-by-feature PCA", "Add features one at a time. Follow how the PC1/PC2 score map changes, without confusing a feature with a principal-component axis.",
          tags$div(class = "pca-step-controls",
            selectInput("pca_order", "Add features in", choices = c(
              "Creation / CSV order" = "created", "Late transition first" = "late",
              "Strongest PC1 loading first" = "loading"), selected = "created"),
            selectInput("pca_mode", "Comparison method", choices = c(
              "Refit and align PCA each step" = "refit",
              "Accumulate contributions in fixed axes" = "fixed"), selected = "refit"),
            sliderInput("pca_step", "Features included", min = 1, max = 12,
                        value = 2, step = 1)),
          tags$div(class = "step-actions",
            actionButton("pca_prev", "− Previous feature", class = "btn-outline-secondary"),
            actionButton("pca_next", "+ Add next feature", class = "btn-outline-primary")),
          uiOutput("pca_step_description"),
          tags$div(class = "lab-grid equal inset-grid",
            plotOutput("pca_progressive_plot", height = "470px"),
            plotOutput("loading_walk_plot", height = "470px")),
          tags$p(class = "hint", "Refit mode aligns each new PCA plane to the full-data plane for visual comparison; the axes are not literally identical. The loading walk always uses the fixed full-data basis."))
      ),
      tags$div(class = "lab-grid one",
        panel("Feature contributions", "PC1 loadings refer to features, not score coordinates.",
          plotOutput("loadings_plot", height = "325px"))
      )
    ),
    bslib::nav_panel("Distance maps",
      tags$div(class = "lab-grid one",
        panel("Distance & embedding workbench",
          "PCA maps measured values; PCoA, NMDS and Sammon map distances. They use the same objects and preprocessing, but only the latter three change when you select a new distance.",
          tags$div(class = "ordination-controls",
            selectInput("distance_method", "Pairwise distance",
              choices = c("Euclidean" = "euclidean", "Manhattan (L1)" = "manhattan",
                          "Canberra" = "canberra", "Bray-Curtis" = "bray",
                          "1 - Pearson correlation" = "correlation"),
              selected = "euclidean"),
            selectInput("distance_orientation", "Objects being mapped",
              choices = c("Time points (across features)" = "time",
                          "Features (across time points)" = "feature"),
              selected = "time"),
            checkboxInput("distance_standardize", "Standardize features (z-score)", value = TRUE),
            sliderInput("distance_cap", "Max sampled time points", min = 40, max = 200,
                        value = 120, step = 20)),
          tags$div(class = "step-actions",
            actionButton("run_ordination", "Recalculate four maps", class = "btn-primary"),
            tags$span(class = "hint", "Bray-Curtis needs nonnegative values: turn off standardization.")),
          uiOutput("ordination_note")
        )
      ),
      tags$div(class = "lab-grid one",
        panel("One-screen comparison", "Four maps fitted to exactly the same selected data. Numeric axes and orientations differ by method; compare relative geometry, not absolute position.",
          plotOutput("ordination_comparison", height = "720px"),
          tableOutput("ordination_quality"),
          uiOutput("ordination_failures"),
          tags$p(class = "hint", "Distance rank = Spearman correlation of pairwise distances; relative error = scale-adjusted Euclidean distance mismatch. NMDS and Sammon also optimize different method-specific stress objectives."))
      ),
      tags$div(class = "lab-grid equal",
        panel("Selected map", "Zoom in on one result from the comparison.",
          selectInput("ordination_focus", "Method", choices = c("PCA", "PCoA", "NMDS", "Sammon")),
          plotOutput("ordination_focus_plot", height = "420px")),
        panel("Distance preservation", "Points near the dashed diagonal retain their relative dissimilarities after a global rescale.",
          plotOutput("ordination_shepard_plot", height = "470px"))
      ),
      tags$div(class = "lab-grid one",
        panel("Input distance matrix", "A readable subset of the exact pairwise distances used for each map.",
          plotOutput("ordination_distance_plot", height = "520px"))
      )
    ),
    bslib::nav_panel("Mean & variance",
      tags$div(class = "lab-grid one",
        panel("Experiment: shift the mean or compress the spread",
          "Perturb one feature inside the selected observation window. This is an illustrative transformation, not a biological treatment simulation.",
          tags$div(class = "pca-step-controls",
            selectInput("perturb_feature", "Feature to perturb", choices = NULL),
            sliderInput("mean_shift", "Requested mean offset", min = -0.30,
                        max = 0.30, value = 0, step = 0.01),
            sliderInput("spread_scale", "Spread multiplier", min = 0.20,
                        max = 2.00, value = 1, step = 0.05)),
          tags$div(class = "experiment-toggle",
            checkboxInput("apply_perturbation", "Apply to PCA, distance maps, correlations and exports", value = FALSE),
            actionButton("reset_perturb", "Reset experiment", class = "btn-outline-secondary")),
          uiOutput("experiment_summary"),
          tags$div(class = "lab-grid equal inset-grid",
            plotOutput("experiment_signal_plot", height = "385px"),
            plotOutput("experiment_spread_plot", height = "385px")),
          tags$p(class = "hint", "The transform is mean(original) + offset + spread × (value − mean(original)), then clipped to [0,1]. Clipping can change both the realized mean and variance."))
      ),
      tags$div(class = "lab-grid equal",
        panel("Before vs after across the time window", "Split the active dataset into earlier and later sampled halves (not independent groups).",
          uiOutput("half_stats"),
          plotOutput("half_change_plot", height = "385px")),
        panel("Selected feature: early vs late", "Mean marker = diamond. Variance reflects within-half time variation.",
          plotOutput("half_spread_plot", height = "425px"))
      )
    ),
    bslib::nav_panel("Correlations",
      tags$div(class = "lab-grid equal",
        panel("Correlation matrix", "Click a tile or choose features on the right.",
          plotOutput("heatmap_plot", height = "500px", click = "heatmap_click")),
        panel("Pair inspector", "Compare a selected pair of features.",
          tags$div(class = "pair-selectors",
            selectInput("feature_a", "Feature A", choices = NULL),
            selectInput("feature_b", "Feature B", choices = NULL)),
          uiOutput("pair_stats"),
          plotOutput("pair_scatter", height = "285px"))
      ),
      tags$div(class = "lab-grid one",
        panel("Pairwise trajectories", "The two signals across your observation window.",
          plotOutput("pair_trajectories", height = "275px"))
      )
    ),
    bslib::nav_panel("Data & export",
      tags$div(class = "lab-grid equal",
        panel("Download reproducible outputs", "Export the underlying data, feature metadata, analyses, and settings.",
          tags$div(class = "download-grid",
            downloadButton("download_data", "Data · CSV"),
            downloadButton("download_metadata", "Metadata · CSV"),
            downloadButton("download_scores", "PCA scores · CSV"),
            downloadButton("download_cor", "Correlations · CSV"),
            downloadButton("download_settings", "Settings · TXT"),
            downloadButton("download_heatmap", "Heatmap · PNG"),
            downloadButton("download_active", "Active window · CSV"),
            downloadButton("download_half_changes", "Mean & variance · CSV"),
            downloadButton("download_step_pca", "PCA step · CSV"),
            downloadButton("download_step_png", "PCA step · PNG"),
            downloadButton("download_ordination_scores", "Distance maps · CSV"),
            downloadButton("download_ordination_quality", "Map quality · CSV"),
            downloadButton("download_ordination_matrix", "Distances · CSV"),
            downloadButton("download_ordination_png", "Map comparison · PNG")
          ),
          tags$p(class = "hint", "Data and metadata exports are original full datasets; analyses and active-window exports reflect the selected window and any enabled perturbation.")),
        panel("Data preview", "First eight rows of the current analysis window.",
          tableOutput("data_preview"))
      ),
      tags$div(class = "lab-grid one",
        panel("Feature provenance", "Generated parameters are saved per feature; imported data have no generator metadata.",
          tableOutput("metadata_preview"))
      )
    ),
    bslib::nav_panel("Guide",
      tags$div(class = "guide-container",
        tags$h2("From time courses to structure"),
        tags$p("Trajectory Lab is an exploratory teaching and prototyping tool. It does not implement EWAS, estimate cell composition, or assign biological significance to its synthetic features."),
        tags$div(class = "guide-grid",
          tags$div(tags$h3("1. Create or import"),
                   tags$p("Generate logistic transitions with reproducible noise, or upload a CSV of bounded beta-scale signals. Settings apply when you press Regenerate dataset.")),
          tags$div(tags$h3("2. Choose a window"),
                   tags$p("All analysis tabs use the same selected time interval. PCA standardizes each non-constant feature within that window.")),
          tags$div(tags$h3("3. Compare features"),
                   tags$p("Pearson or rank-based Spearman correlation is computed across time points. Optional Gaussian weights emphasize a chosen time region.")),
          tags$div(tags$h3("4. Export and verify"),
                   tags$p("Download datasets, mean–variance comparisons, metadata, scores, correlations, and settings. Use the supplied tests before deploying changes."))
        ),
        tags$div(class = "guide-grid",
          tags$div(tags$h3("5. Add features step by step"),
                   tags$p("The PCA stepper adds measurement features, not extra PC axes. Adding PC3 to an unchanged PCA cannot move its PC1/PC2 scores. Refit and align each step, or accumulate contributions in a fixed PCA basis. The loading walk is separate from observation scores.")),
          tags$div(tags$h3("6. Change mean and spread"),
                   tags$p("A controlled shift/scale transformation lets you distinguish location from dispersion. The early–late comparison is descriptive, not a two-group hypothesis test.")),
          tags$div(tags$h3("7. Compare embeddings"),
                   tags$p("PCA uses the measurement matrix; PCoA, NMDS, and Sammon use the dissimilarity matrix. Try Euclidean, Manhattan, Canberra, Bray-Curtis, or correlation distance.")),
          tags$div(tags$h3("8. Read distance diagnostics"),
                   tags$p("Pairwise distance rank and scale-adjusted error are shared diagnostics. PCoA can have negative eigenvalues; NMDS and Sammon require positive off-diagonal distances."))
        ),
        tags$div(class = "guide-note",
          tags$strong("Interpretation matters."),
          tags$p("Clipping the simulator's beta-like values to 0–1 and scaling for PCA may influence the result. Unequal time sampling, dependence between observations, latent confounding, and measurement design are not addressed by this demo."))
      )
    )
  ),
  tags$footer(class = "site-footer",
              "TRAJECTORY LAB · Open-source Shiny prototype · Synthetic signal exploration")
)

server <- function(input, output, session) {
  simulated_state <- reactiveVal(simulate_dataset())

  current_dataset <- reactive({
    if (identical(input$source, "uploaded")) {
      req(input$user_csv)
      result <- safe_result(read_beta_csv(input$user_csv$datapath))
      validate(need(!inherits(result, "lab_error"),
                    if (inherits(result, "lab_error")) result$message else TRUE))
      result
    } else {
      simulated_state()
    }
  })

  observeEvent(input$generate, {
    new_data <- safe_result(simulate_dataset(
      n_features = input$n_features, n_timepoints = input$n_timepoints,
      pattern = input$pattern, amplitude = input$amplitude,
      transition_days = input$transition_days, noise_sd = input$noise_sd,
      seed = input$seed
    ))
    if (inherits(new_data, "lab_error")) {
      showNotification(new_data$message, type = "error")
    } else {
      simulated_state(new_data)
      showNotification("Dataset generated with the specified seed.", type = "message", duration = 3)
    }
  })

  observeEvent(input$add_feature, {
    updated <- safe_result(add_custom_feature(
      simulated_state(), midpoint = input$custom_midpoint,
      signed_amplitude = input$custom_amplitude, baseline = input$custom_baseline,
      transition_days = input$custom_width,
      noise_sd = input$noise_sd
    ))
    if (inherits(updated, "lab_error")) {
      showNotification(updated$message, type = "warning")
    } else {
      simulated_state(updated)
      updateSelectInput(session, "selected_feature", selected = tail(feature_names(updated$data), 1L))
    }
  })

  observeEvent(input$undo_feature, {
    updated <- safe_result(undo_custom_feature(simulated_state()))
    if (inherits(updated, "lab_error")) {
      showNotification(updated$message, type = "warning")
    } else {
      simulated_state(updated)
    }
  })

  observeEvent(list(input$source, input$user_csv), {
    if (identical(input$source, "uploaded") && !is.null(input$user_csv)) {
      data <- safe_result(read_beta_csv(input$user_csv$datapath))
      if (inherits(data, "lab_error")) {
        showNotification(data$message, type = "error", duration = 8)
      } else {
        bounds <- range(data$data$time)
        updateSliderInput(session, "time_window", min = bounds[1L], max = bounds[2L],
                          value = bounds, step = diff(bounds) / 200)
        updateNumericInput(session, "focus_time", value = mean(bounds))
        updateNumericInput(session, "bandwidth", value = max(0.000001, diff(bounds) / 6))
      }
    } else if (identical(input$source, "simulated")) {
      updateSliderInput(session, "time_window", min = 0, max = 5000,
                        value = c(0, 5000), step = 50)
      updateNumericInput(session, "focus_time", value = 2500)
      updateNumericInput(session, "bandwidth", value = 750)
    }
  }, ignoreInit = TRUE)

  observeEvent(current_dataset(), {
    features <- feature_names(current_dataset()$data)
    choose <- function(old, default) if (length(old) == 1L && old %in% features) old else default
    updateSelectInput(session, "selected_feature", choices = features,
                      selected = choose(isolate(input$selected_feature), features[1L]))
    updateSelectInput(session, "feature_a", choices = features,
                      selected = choose(isolate(input$feature_a), features[1L]))
    updateSelectInput(session, "feature_b", choices = features,
                      selected = choose(isolate(input$feature_b), features[min(2L, length(features))]))
    updateSelectInput(session, "perturb_feature", choices = features,
                      selected = choose(isolate(input$perturb_feature), features[1L]))
  })

  windowed_data <- reactive({
    data <- safe_result(window_data(current_dataset()$data, input$time_window))
    validate(need(!inherits(data, "lab_error"),
                  if (inherits(data, "lab_error")) data$message else TRUE))
    data
  })

  previous_feature_count <- reactiveVal(NULL)
  observeEvent(windowed_data(), {
    total <- length(active_features(windowed_data()))
    if (total >= 2L) {
      requested <- isolate(input$pca_step)
      last <- isolate(previous_feature_count())
      if (is.null(last) || total > last) requested <- total
      if (length(requested) != 1L || is.na(requested)) requested <- 2L
      minimum <- if (identical(isolate(input$pca_mode), "fixed")) 1L else 2L
      updateSliderInput(session, "pca_step", min = minimum, max = total,
                        value = min(total, max(minimum, as.integer(requested))))
    }
    previous_feature_count(total)
  })

  observeEvent(input$pca_mode, {
    total <- length(active_features(windowed_data()))
    if (total >= 2L) {
      minimum <- if (identical(input$pca_mode, "fixed")) 1L else 2L
      requested <- isolate(input$pca_step)
      if (length(requested) != 1L || is.na(requested)) requested <- minimum
      updateSliderInput(session, "pca_step", min = minimum, max = total,
                        value = min(total, max(minimum, as.integer(requested))))
    }
  })

  observeEvent(input$reset_perturb, {
    updateSliderInput(session, "mean_shift", value = 0)
    updateSliderInput(session, "spread_scale", value = 1)
    updateCheckboxInput(session, "apply_perturbation", value = FALSE)
  })

  experiment_preview <- reactive({
    req(input$perturb_feature)
    result <- safe_result(perturb_feature(
      windowed_data(), input$perturb_feature,
      shift = input$mean_shift, spread = input$spread_scale))
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  })

  analyzed_data <- reactive({
    data <- windowed_data()
    if (isTRUE(input$apply_perturbation) &&
        length(input$perturb_feature) == 1L &&
        input$perturb_feature %in% feature_names(data)) {
      data <- experiment_preview()$data
    }
    data
  })

  pca_result <- reactive(safe_result(run_pca(analyzed_data())))
  cor_result <- reactive(safe_result(run_correlation(
    analyzed_data(), method = input$cor_method,
    weighted = isTRUE(input$weighted),
    focus = input$focus_time, bandwidth = input$bandwidth
  )))

  checked_pca <- reactive({
    result <- pca_result()
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  })
  checked_cor <- reactive({
    result <- cor_result()
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  })

  output$metrics <- renderUI({
    data <- analyzed_data()
    pca <- pca_result()
    cor <- cor_result()
    pca_share <- if (inherits(pca, "lab_error")) "N/A" else sprintf("%.1f%%", 100 * sum(head(pca$explained$fraction, 2L)))
    cor_mean <- if (inherits(cor, "lab_error") || !is.finite(cor$mean_abs)) "N/A" else sprintf("%.3f", cor$mean_abs)
    tags$div(class = "metric-grid",
      metric_tile("FEATURES", length(feature_names(data)),
                  if (identical(current_dataset()$source, "simulated")) "Synthetic profiles" else "Uploaded profiles"),
      metric_tile("TIME POINTS", format(nrow(data), big.mark = ","),
                  if (isTRUE(input$apply_perturbation)) "Perturbation active" else "Within selected window"),
      metric_tile("PC1 + PC2", pca_share,
                  "Explained variance"),
      metric_tile("MEAN |r|", cor_mean,
                  "Unique off-diagonal pairs")
    )
  })

  output$trajectory_plot <- renderPlot({
    data <- analyzed_data()
    req(input$selected_feature)
    validate(need(input$selected_feature %in% feature_names(data), "Select a feature."))
    plot_trajectories(data, input$selected_feature)
  }, res = 110)

  output$feature_stats <- renderUI({
    data <- current_dataset()
    req(input$selected_feature)
    meta <- data$metadata[data$metadata$feature == input$selected_feature, , drop = FALSE]
    validate(need(nrow(meta) > 0L, "Select a feature."))
    if (is.finite(meta$midpoint[1L])) {
      tags$div(class = "feature-statline",
        tags$span(sprintf("Midpoint: %s", format(round(meta$midpoint[1L]), big.mark = ","))),
        tags$span(sprintf("Signed amplitude: %+0.2f", meta$signed_amplitude[1L]))
      )
    } else {
      tags$div(class = "feature-statline", tags$span("Imported feature · no generator parameters"))
    }
  })

  output$feature_plot <- renderPlot({
    data <- analyzed_data()
    req(input$selected_feature)
    validate(need(input$selected_feature %in% feature_names(data), "Select a feature."))
    meta <- current_dataset()$metadata
    plot_feature(data, input$selected_feature,
                 meta[meta$feature == input$selected_feature, , drop = FALSE])
  }, res = 110)

  output$pca_plot <- renderPlot(plot_pca_scores(checked_pca()), res = 110)
  output$scree_plot <- renderPlot(plot_pca_scree(checked_pca()), res = 110)
  output$loadings_plot <- renderPlot(plot_pca_loadings(checked_pca()), res = 110)

  pca_sequence <- reactive({
    req(input$pca_order)
    result <- safe_result(feature_sequence(
      analyzed_data(), current_dataset()$metadata, input$pca_order))
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  })
  pca_step_result <- reactive({
    req(input$pca_mode, input$pca_step)
    sequence <- pca_sequence()
    min_step <- if (identical(input$pca_mode, "fixed")) 1L else 2L
    step <- min(length(sequence), max(min_step, as.integer(input$pca_step)))
    result <- safe_result(progressive_pca(
      analyzed_data(), sequence = sequence, step = step, mode = input$pca_mode))
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  })
  observeEvent(input$pca_prev, {
    total <- length(pca_sequence())
    old <- isolate(input$pca_step)
    if (length(old) == 1L && is.finite(old))
      updateSliderInput(session, "pca_step", value = max(
        if (identical(isolate(input$pca_mode), "fixed")) 1L else 2L,
        min(total, old - 1L)))
  })
  observeEvent(input$pca_next, {
    total <- length(pca_sequence())
    old <- isolate(input$pca_step)
    if (length(old) == 1L && is.finite(old))
      updateSliderInput(session, "pca_step", value = min(total, max(1L, old + 1L)))
  })
  output$pca_step_description <- renderUI({
    p <- pca_step_result()
    tags$div(class = "feature-statline",
      tags$span(sprintf("Step %d/%d", p$step, length(p$sequence))),
      tags$span(paste("Just added:", tail(p$included, 1L))),
      tags$span(if (identical(p$mode, "refit"))
                  "New PCA refit + orthogonal alignment" else
                  "Full PCA basis stays fixed"))
  })
  output$pca_progressive_plot <- renderPlot(plot_progressive_pca(pca_step_result()), res = 105)
  output$loading_walk_plot <- renderPlot({
    p <- pca_step_result()
    plot_loading_walk(loading_walk(analyzed_data(), p$sequence, p$step))
  }, res = 105)

  half_summary <- reactive(half_window_summary(analyzed_data()))
  output$half_stats <- renderUI({
    req(input$perturb_feature)
    summary <- half_summary()
    row <- summary[summary$feature == input$perturb_feature, , drop = FALSE]
    validate(need(nrow(row) == 1L, "Select a feature."))
    ratio <- if (is.finite(row$variance_ratio))
      sprintf("%.3fx", row$variance_ratio) else "N/A (zero baseline variance)"
    tags$div(class = "feature-statline",
      tags$span(sprintf("Early mean %.3f", row$mean_early)),
      tags$span(sprintf("Late mean %.3f", row$mean_late)),
      tags$span(sprintf("Δ mean %+.3f", row$delta_mean)),
      tags$span(paste("Variance ratio", ratio)))
  })
  output$half_change_plot <- renderPlot({
    req(input$perturb_feature)
    result <- safe_result(plot_window_change_map(half_summary(), input$perturb_feature))
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  }, res = 105)
  output$half_spread_plot <- renderPlot({
    req(input$perturb_feature)
    plot_half_distributions(analyzed_data(), input$perturb_feature)
  }, res = 105)
  output$experiment_summary <- renderUI({
    e <- experiment_preview()
    variance <- if (is.finite(e$variance_ratio))
      sprintf("%.3fx", e$variance_ratio) else "N/A (zero baseline variance)"
    tags$div(class = "feature-statline",
      tags$span(sprintf("Actual mean shift: %+.3f", e$actual_shift)),
      tags$span(paste("Variance ratio:", variance)),
      tags$span(sprintf("Clipped: %d/%d", e$clipped, e$total)),
      if (isTRUE(input$apply_perturbation))
        tags$span("Experiment ACTIVE in analysis") else
        tags$span("Preview only"))
  })
  output$experiment_signal_plot <- renderPlot({
    plot_perturbed_trajectory(experiment_preview()$paired, input$perturb_feature)
  }, res = 105)
  output$experiment_spread_plot <- renderPlot({
    plot_perturbed_spread(experiment_preview()$paired)
  }, res = 105)

  ordination_result <- eventReactive(input$run_ordination, {
    result <- safe_result(run_ordination_suite(
      analyzed_data(),
      distance = input$distance_method,
      standardize = isTRUE(input$distance_standardize),
      orientation = input$distance_orientation,
      max_points = as.integer(input$distance_cap)
    ))
    validate(need(!inherits(result, "lab_error"),
                  if (inherits(result, "lab_error")) result$message else TRUE))
    result
  }, ignoreNULL = FALSE)

  checked_ordination <- reactive(ordination_result())

  output$ordination_note <- renderUI({
    suite <- checked_ordination()
    tags$div(class = "feature-statline",
      tags$span(sprintf("Objects in maps: %d of %d", suite$data$n_objects,
                        suite$data$total_objects)),
      tags$span(sprintf("Measurements per object: %d", suite$data$n_measurements)),
      tags$span(paste("Distance:", distance_labels[[suite$distance_method]])),
      if (length(suite$data$features_removed))
        tags$span(sprintf("Constant features omitted: %d", length(suite$data$features_removed))))
  })
  output$ordination_comparison <- renderPlot(
    plot_ordination_comparison(checked_ordination()), res = 105)
  output$ordination_quality <- renderTable({
    d <- checked_ordination()$summary
    data.frame(Method = d$method,
               `Distance rank` = ifelse(is.finite(d$distance_rank),
                                        sprintf("%.3f", d$distance_rank), "N/A"),
               `Relative error` = ifelse(is.finite(d$distance_error),
                                         sprintf("%.3f", d$distance_error), "N/A"),
               Status = d$status, check.names = FALSE)
  }, striped = TRUE, hover = TRUE, spacing = "xs")
  output$ordination_failures <- renderUI({
    d <- checked_ordination()$summary
    unavailable <- d[d$status != "OK", , drop = FALSE]
    if (!nrow(unavailable)) return(NULL)
    tags$div(class = "hint", tags$strong("Unavailable maps: "),
      paste(paste(unavailable$method, unavailable$detail, sep = ": "), collapse = " | "))
  })
  valid_ordination_focus <- reactive({
    req(input$ordination_focus)
    suite <- checked_ordination()
    fit <- suite$fits[[input$ordination_focus]]
    validate(need(!is.null(fit) && is.null(fit$error),
                  if (is.null(fit)) "Choose an available method." else
                    if (!is.null(fit$error)) fit$error else TRUE))
    input$ordination_focus
  })
  output$ordination_focus_plot <- renderPlot({
    plot_ordination_focus(checked_ordination(), valid_ordination_focus())
  }, res = 105)
  output$ordination_shepard_plot <- renderPlot({
    plot_ordination_shepard(checked_ordination(), valid_ordination_focus())
  }, res = 105)
  output$ordination_distance_plot <- renderPlot({
    plot_dissimilarity_matrix(checked_ordination())
  }, res = 105)

  output$heatmap_plot <- renderPlot(plot_correlation_heatmap(checked_cor()), res = 110)

  observeEvent(input$heatmap_click, {
    features <- feature_names(analyzed_data())
    click <- input$heatmap_click
    x <- as.integer(round(click$x))
    y <- as.integer(round(click$y))
    if (length(x) == 1L && length(y) == 1L &&
        is.finite(x) && is.finite(y) &&
        x >= 1L && x <= length(features) && y >= 1L && y <= length(features)) {
      updateSelectInput(session, "feature_a", selected = features[x])
      updateSelectInput(session, "feature_b", selected = rev(features)[y])
    }
  })

  selected_pair <- reactive({
    features <- feature_names(analyzed_data())
    req(input$feature_a, input$feature_b)
    validate(need(input$feature_a %in% features && input$feature_b %in% features,
                  "Select two features."))
    c(input$feature_a, input$feature_b)
  })

  output$pair_stats <- renderUI({
    pair <- selected_pair()
    cor <- checked_cor()
    value <- cor$matrix[pair[1L], pair[2L]]
    tags$div(class = "feature-statline",
      tags$span(if (is.finite(value)) sprintf("Correlation: %.3f", value) else "Correlation: N/A"),
      tags$span(sprintf("Effective N: %.1f", cor$effective_n))
    )
  })
  output$pair_scatter <- renderPlot({
    pair <- selected_pair()
    plot_pair_scatter(analyzed_data(), pair[1L], pair[2L])
  }, res = 110)
  output$pair_trajectories <- renderPlot({
    pair <- selected_pair()
    plot_pair_trajectories(analyzed_data(), pair[1L], pair[2L])
  }, res = 110)

  output$data_preview <- renderTable({
    data <- head(analyzed_data(), 8L)
    data[, head(names(data), 9L), drop = FALSE]
  }, striped = TRUE, hover = TRUE, digits = 3, spacing = "xs")

  output$metadata_preview <- renderTable({
    head(current_dataset()$metadata[, c("feature", "midpoint", "signed_amplitude",
                                        "transition_days", "source")], 40L)
  }, striped = TRUE, hover = TRUE, digits = 3, spacing = "xs")

  output$download_data <- downloadHandler(
    filename = function() "trajectory-lab-data.csv",
    content = function(file) utils::write.csv(current_dataset()$data, file, row.names = FALSE)
  )
  output$download_metadata <- downloadHandler(
    filename = function() "trajectory-lab-metadata.csv",
    content = function(file) utils::write.csv(current_dataset()$metadata, file, row.names = FALSE,
                                                na = "")
  )
  output$download_scores <- downloadHandler(
    filename = function() "trajectory-lab-pca-scores.csv",
    content = function(file) utils::write.csv(checked_pca()$scores, file, row.names = FALSE)
  )
  output$download_cor <- downloadHandler(
    filename = function() "trajectory-lab-correlation.csv",
    content = function(file) {
      mat <- checked_cor()$matrix
      utils::write.csv(data.frame(feature = rownames(mat), mat, check.names = FALSE),
                       file, row.names = FALSE, na = "")
    }
  )
  output$download_settings <- downloadHandler(
    filename = function() "trajectory-lab-settings.txt",
    content = function(file) writeLines(c(make_session_notes(
      current_dataset(), input$time_window, input$cor_method,
      isTRUE(input$weighted), input$focus_time, input$bandwidth),
      "", "PCA STEP-THROUGH:",
      paste("Ordering:", input$pca_order),
      paste("Step:", input$pca_step),
      paste("Mode:", input$pca_mode),
      "", "SHIFT / SPREAD EXPERIMENT:",
      paste("Feature:", input$perturb_feature),
      paste("Requested mean offset:", input$mean_shift),
      paste("Spread multiplier:", input$spread_scale),
      paste("Applied to analysis:", isTRUE(input$apply_perturbation)),
      "", "DISTANCE MAPS (most recently calculated):",
      paste("Distance metric:", checked_ordination()$distance_method),
      paste("Objects:", checked_ordination()$data$orientation),
      paste("Standardized:", checked_ordination()$data$scaled),
      paste("Mapped objects:", checked_ordination()$data$n_objects),
      paste("Max sampled time points:", checked_ordination()$data$max_points)
    ), file)
  )
  output$download_active <- downloadHandler(
    filename = function() "trajectory-lab-active-window.csv",
    content = function(file) utils::write.csv(analyzed_data(), file, row.names = FALSE)
  )
  output$download_half_changes <- downloadHandler(
    filename = function() "trajectory-lab-mean-variance.csv",
    content = function(file) utils::write.csv(half_summary(), file, row.names = FALSE)
  )
  output$download_step_pca <- downloadHandler(
    filename = function() "trajectory-lab-pca-step.csv",
    content = function(file) {
      p <- pca_step_result()
      d <- data.frame(time = p$current$time,
                      current_PC1 = p$current$PC1, current_PC2 = p$current$PC2,
                      reference_PC1 = p$reference$PC1,
                      reference_PC2 = p$reference$PC2)
      utils::write.csv(d, file, row.names = FALSE)
    }
  )
  output$download_step_png <- downloadHandler(
    filename = function() "trajectory-lab-pca-step.png",
    content = function(file) ggplot2::ggsave(
      filename = file, plot = plot_progressive_pca(pca_step_result()),
      width = 9, height = 7, dpi = 160, units = "in", device = "png")
  )
  output$download_ordination_scores <- downloadHandler(
    filename = function() "trajectory-lab-ordination-scores.csv",
    content = function(file) utils::write.csv(ordination_scores(checked_ordination()),
                                              file, row.names = FALSE)
  )
  output$download_ordination_quality <- downloadHandler(
    filename = function() "trajectory-lab-ordination-quality.csv",
    content = function(file) utils::write.csv(checked_ordination()$summary,
                                              file, row.names = FALSE)
  )
  output$download_ordination_matrix <- downloadHandler(
    filename = function() "trajectory-lab-distance-matrix.csv",
    content = function(file) {
      mat <- as.matrix(checked_ordination()$distance)
      utils::write.csv(data.frame(object = rownames(mat), mat, check.names = FALSE),
                       file, row.names = FALSE)
    }
  )
  output$download_ordination_png <- downloadHandler(
    filename = function() "trajectory-lab-map-comparison.png",
    content = function(file) ggplot2::ggsave(
      filename = file, plot = plot_ordination_comparison(checked_ordination()),
      width = 11, height = 9, dpi = 160, units = "in", device = "png")
  )
  output$download_heatmap <- downloadHandler(
    filename = function() "trajectory-lab-heatmap.png",
    content = function(file) ggplot2::ggsave(
      filename = file, plot = plot_correlation_heatmap(checked_cor()),
      width = 8, height = 7, dpi = 160, units = "in", device = "png")
  )
}

shinyApp(ui = ui, server = server)
