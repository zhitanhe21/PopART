#!/usr/bin/env Rscript

# Regenerate the repository-only synthetic teaching CSV files. Run from the
# repository root. This script is simulation tooling and is not installed.

project_dir <- normalizePath(getwd(), mustWork = TRUE)
loader <- file.path(project_dir, "simulation", "R", "load_simulation.R")
if (!file.exists(loader) || !file.exists(file.path(project_dir, "DESCRIPTION"))) {
  stop("Run this script from the popart repository root.", call. = FALSE)
}

source(loader)
study <- new.env(parent = globalenv())
source_reference_study(project_dir, envir = study)
generated <- study$generate_reference_data(
  m = 20L,
  n_trial = 200L,
  n_auxiliary = 200L,
  p_resp = 0.5,
  p_cens = 0.3,
  seed = 20260805L
)
data <- generated$dat

trial <- data[
  data$S == 1L,
  c("id", "clust", "X1", "W1", "W2", "W3", "A", "R", "C", "Y")
]
names(trial) <- c(
  "participant_id",
  "cluster",
  "baseline_risk",
  "baseline_binary",
  "baseline_score_1",
  "baseline_score_2",
  "treatment",
  "responded",
  "censored",
  "outcome"
)
trial$baseline_score_1 <- round(trial$baseline_score_1, 6L)
trial$baseline_score_2 <- round(trial$baseline_score_2, 6L)
trial$censored[trial$responded == 0L] <- NA_integer_
trial$outcome[trial$responded == 0L | trial$censored == 1L] <- NA_integer_

auxiliary <- data[
  data$S == 0L,
  c("id", "clust", "X1", "W1", "W2", "W3", "wt")
]
names(auxiliary) <- c(
  "participant_id",
  "cluster",
  "baseline_risk",
  "baseline_binary",
  "baseline_score_1",
  "baseline_score_2",
  "survey_weight"
)
auxiliary$baseline_score_1 <- round(auxiliary$baseline_score_1, 6L)
auxiliary$baseline_score_2 <- round(auxiliary$baseline_score_2, 6L)
auxiliary$survey_weight <- round(auxiliary$survey_weight, 6L)

output_dir <- file.path(project_dir, "simulation", "data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  trial,
  file.path(output_dir, "example_trial.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  auxiliary,
  file.path(output_dir, "example_auxiliary.csv"),
  row.names = FALSE,
  na = ""
)
cat("Regenerated repository teaching data in", output_dir, "\n")
