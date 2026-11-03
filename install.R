# Run once in R from the repository root: source("install.R")
required <- c("shiny", "bslib", "ggplot2", "MASS")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
message("Dependencies available: ", paste(required, collapse = ", "))
