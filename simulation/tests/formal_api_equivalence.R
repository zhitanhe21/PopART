#!/usr/bin/env Rscript

# Verify that the formal two-data-frame API reproduces the preserved reference
# equations and HAL fold assignments for one deterministic study replicate.

find_repository_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(current, "simulation", "R", "load_simulation.R"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the popart repository root.", call. = FALSE)
    }
    current <- parent
  }
}


if (!requireNamespace("popart", quietly = TRUE)) {
  stop("Install popart before running this equivalence test.", call. = FALSE)
}

project_dir <- find_repository_root()
source(file.path(project_dir, "simulation", "R", "load_simulation.R"))
study <- new.env(parent = globalenv())
source_reference_study(project_dir, envir = study)

seed <- 91000001L
n_trial <- 120L
n_auxiliary <- 120L
generated <- study$generate_reference_data(
  m = 20L,
  n_trial = n_trial,
  n_auxiliary = n_auxiliary,
  p_resp = 0.5,
  p_cens = 0.3,
  seed = seed
)

reference_config <- list(
  cv_folds = 3L,
  cv_nlambda = 3L,
  cv_workers_per_fit = 1L,
  fit_workers = 1L,
  backend = "auto"
)
reference <- suppressWarnings(study$run_reference_replicate(
  m = 20L,
  n_trial = n_trial,
  n_aux = n_auxiliary,
  p_resp = 0.5,
  p_cens = 0.3,
  seed = seed,
  config = reference_config,
  use_cache = FALSE
))

formal_control <- popart::popart_control(
  # The preserved reference estimator uses full-sample nuisance predictions.
  # Disable the package's default outer cross-fitting for this equivalence test.
  crossfit_folds = 1L,
  n_cv_folds = 3L,
  n_lambda_values = 3L,
  n_cv_workers = 1L,
  n_fit_workers = 1L,
  random_seed = seed,
  num_knots = c(5L, 3L),
  normalize_auxiliary_weights = FALSE
)
formal <- suppressWarnings(study$fit_generated_data_with_popart(
  generated,
  control = formal_control
))

reference_rows <- reference$results
reference_estimates <- c(
  reference_rows$etahat_0[1L],
  reference_rows$etahat_1[1L],
  reference_rows$rdhat[1L],
  reference_rows$rrhat[1L],
  reference_rows$etahat_0[2L],
  reference_rows$etahat_1[2L],
  reference_rows$rdhat[2L],
  reference_rows$rrhat[2L]
)
reference_variances <- c(
  reference_rows$cov_00[1L],
  reference_rows$cov_11[1L],
  reference_rows$var_rd[1L],
  reference_rows$var_rr[1L],
  reference_rows$cov_00[2L],
  reference_rows$cov_11[2L],
  reference_rows$var_rd[2L],
  reference_rows$var_rr[2L]
)

estimate_difference <- max(abs(
  reference_estimates - formal$estimates$estimate
))
variance_difference <- max(abs(
  reference_variances - formal$variance$variance
))
stopifnot(
  isTRUE(all.equal(
    reference_estimates,
    formal$estimates$estimate,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    reference_variances,
    formal$variance$variance,
    tolerance = 1e-12
  ))
)

cat("Formal API/reference equivalence passed.\n")
cat("  maximum estimate difference:", format(estimate_difference), "\n")
cat("  maximum variance difference:", format(variance_difference), "\n")
