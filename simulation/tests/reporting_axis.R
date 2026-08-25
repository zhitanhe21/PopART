#!/usr/bin/env Rscript

# Repository-only smoke test for the manuscript-style nested sample-size axis.

find_repository_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    marker <- file.path(current, "simulation", "R", "load_simulation.R")
    if (file.exists(marker) && file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Could not locate the popart repository root.", call. = FALSE)
}

project_dir <- find_repository_root()
source(file.path(project_dir, "simulation", "R", "load_simulation.R"))
study <- new.env(parent = globalenv())
source_reference_study(project_dir, envir = study)
study$load_reference_study_packages()

expected <- c("200.200", "500.200", "200.500", "500.500")
axis <- study$make_sample_size_axis(
  n_aux = c(200, 500, 200, 500),
  n_trial = c(200, 200, 500, 500)
)
stopifnot(identical(levels(axis), expected))

set.seed(17)
order <- sample.int(length(axis))
shuffled <- study$make_sample_size_axis(
  n_aux = c(200, 500, 200, 500)[order],
  n_trial = c(200, 200, 500, 500)[order]
)
stopifnot(identical(levels(shuffled), expected))

grid <- expand.grid(
  seed = 1:3,
  Version = c("Naive", "Proposed"),
  n_trial = c(200L, 500L),
  n_aux = c(200L, 500L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
grid$Estimator <- "AIPW"
grid$eta_0 <- 0.4
grid$eta_1 <- 0.3
grid$rd <- -0.1
grid$rr <- 0.75
offset <- (grid$seed - 2) * 0.01
grid$etahat_0 <- grid$eta_0 + offset
grid$etahat_1 <- grid$eta_1 - offset
grid$rdhat <- grid$etahat_1 - grid$etahat_0
grid$rrhat <- grid$etahat_1 / grid$etahat_0
grid$cov_00 <- 0.002
grid$cov_01 <- 0.0005
grid$cov_11 <- 0.002
grid$var_rd <- grid$cov_11 + grid$cov_00 - 2 * grid$cov_01
grid$var_rr <- (
  grid$cov_11 / grid$etahat_0^2 +
    grid$cov_00 * grid$etahat_1^2 / grid$etahat_0^4 -
    2 * grid$cov_01 * grid$etahat_1 / grid$etahat_0^3
)

run_id <- "reporting_axis_smoke"
out_dir <- file.path(tempdir(), run_id)
figure_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  grid,
  file.path(out_dir, paste0("monte_carlo_", run_id, "_results.csv")),
  row.names = FALSE
)

paths <- study$write_monte_carlo_figures(
  run_id = run_id,
  out_dir = out_dir,
  figure_dir = figure_dir,
  width = 4,
  height = 3,
  dpi = 72
)
stopifnot(
  identical(names(paths), c("estimates", "variance", "confidence")),
  all(file.exists(paths)),
  all(file.info(paths)$size > 0)
)

cat("Simulation reporting-axis smoke test passed.\n")
