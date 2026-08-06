###############################################################################
# User-facing estimation controls
###############################################################################


#' Configure PopART estimation
#'
#' Creates a validated control object for the HAL nuisance fits and Wald
#' inference used by [fit_popart()]. The defaults preserve the current
#' computational settings while keeping all resource choices explicit.
#'
#' @param n_cv_folds Number of cross-validation folds used by each HAL fit.
#' @param n_lambda_values Number of lasso penalty values evaluated within each fit.
#' @param n_cv_workers Number of workers used across cross-validation folds when
#'   complete nuisance fits are run serially.
#' @param n_fit_workers Number of complete nuisance fits run concurrently. Values
#'   greater than one require `n_cv_workers = 1` to avoid nested parallelism.
#' @param parallel_backend Complete-fit parallel backend: `"auto"`, `"fork"`, or
#'   `"psock"`.
#' @param random_seed Nonnegative integer used to construct reproducible CV folds.
#' @param smoothness_orders,max_degree,num_knots HAL basis settings.
#' @param treatment_probability Known probability of assignment to the treated
#'   arm. The default corresponds to equal randomization.
#' @param conf_level Confidence level for Wald intervals.
#' @param normalize_auxiliary_weights Whether auxiliary weights should be
#'   rescaled to have mean one. This is the scale used by the reference
#'   estimator and makes their scale relative to unit-weighted trial records
#'   explicit.
#' @param keep_nuisance_fits Whether fitted HAL objects are retained in the
#'   returned object.
#' @param positivity_threshold Threshold used only to flag small estimated
#'   probabilities in diagnostics. Predictions are not truncated.
#'
#' @return An object of class `popart_control`.
#' @export
popart_control <- function(
    n_cv_folds = 3L,
    n_lambda_values = 50L,
    n_cv_workers = 1L,
    n_fit_workers = 1L,
    parallel_backend = c("auto", "fork", "psock"),
    random_seed = 1L,
    smoothness_orders = 1L,
    max_degree = 2L,
    num_knots = c(5L, 3L),
    treatment_probability = 0.5,
    conf_level = 0.95,
    normalize_auxiliary_weights = TRUE,
    keep_nuisance_fits = FALSE,
    positivity_threshold = 0.01) {
  parallel_backend <- match.arg(parallel_backend)

  n_cv_folds <- .popart_count(n_cv_folds, "n_cv_folds", minimum = 3L)
  n_lambda_values <- .popart_count(n_lambda_values, "n_lambda_values")
  n_cv_workers <- .popart_count(n_cv_workers, "n_cv_workers")
  n_fit_workers <- .popart_count(
    n_fit_workers,
    "n_fit_workers",
    maximum = 4L
  )
  random_seed <- .popart_count(
    random_seed,
    "random_seed",
    minimum = 0L,
    maximum = .Machine$integer.max - 100000L
  )
  max_degree <- .popart_count(max_degree, "max_degree")

  if (n_cv_workers > n_cv_folds) {
    stop("n_cv_workers cannot exceed n_cv_folds.", call. = FALSE)
  }
  if (n_fit_workers > 1L && n_cv_workers > 1L) {
    stop(
      "Use either parallel complete fits (n_fit_workers > 1) or parallel CV ",
      "(n_cv_workers > 1), not both.",
      call. = FALSE
    )
  }
  if (identical(parallel_backend, "fork") && .Platform$OS.type == "windows") {
    stop("The fork backend is not available on Windows.", call. = FALSE)
  }

  smoothness_orders <- .popart_integer_vector(
    smoothness_orders,
    "smoothness_orders",
    minimum = 0L
  )
  num_knots <- .popart_integer_vector(num_knots, "num_knots", minimum = 1L)

  treatment_probability <- .popart_probability(
    treatment_probability,
    "treatment_probability"
  )
  conf_level <- .popart_probability(
    conf_level,
    "conf_level"
  )
  positivity_threshold <- .popart_probability(
    positivity_threshold,
    "positivity_threshold"
  )
  if (positivity_threshold >= 0.5) {
    stop("positivity_threshold must be less than 0.5.", call. = FALSE)
  }

  normalize_auxiliary_weights <- .popart_flag(
    normalize_auxiliary_weights,
    "normalize_auxiliary_weights"
  )
  keep_nuisance_fits <- .popart_flag(
    keep_nuisance_fits,
    "keep_nuisance_fits"
  )

  structure(
    list(
      n_cv_folds = n_cv_folds,
      n_lambda_values = n_lambda_values,
      n_cv_workers = n_cv_workers,
      n_fit_workers = n_fit_workers,
      parallel_backend = parallel_backend,
      random_seed = random_seed,
      smoothness_orders = smoothness_orders,
      max_degree = max_degree,
      num_knots = num_knots,
      treatment_probability = treatment_probability,
      conf_level = conf_level,
      normalize_auxiliary_weights = normalize_auxiliary_weights,
      keep_nuisance_fits = keep_nuisance_fits,
      positivity_threshold = positivity_threshold
    ),
    class = "popart_control"
  )
}


#' @export
print.popart_control <- function(x, ...) {
  cat("PopART estimation controls\n")
  cat("  CV folds / workers: ", x$n_cv_folds, " / ", x$n_cv_workers, "\n", sep = "")
  cat("  Complete-fit workers: ", x$n_fit_workers, "\n", sep = "")
  cat("  Parallel backend: ", x$parallel_backend, "\n", sep = "")
  cat("  HAL knots: ", paste(x$num_knots, collapse = ", "), "\n", sep = "")
  cat("  Confidence level: ", x$conf_level, "\n", sep = "")
  invisible(x)
}


.popart_count <- function(x, name, minimum = 1L,
                          maximum = .Machine$integer.max) {
  value <- suppressWarnings(as.integer(x))
  if (length(x) != 1L || is.na(value) || !is.finite(as.numeric(x)) ||
      as.numeric(x) != value || value < minimum || value > maximum) {
    stop(
      name, " must be an integer between ", minimum, " and ", maximum, ".",
      call. = FALSE
    )
  }
  value
}


.popart_integer_vector <- function(x, name, minimum = 0L) {
  value <- suppressWarnings(as.integer(x))
  if (length(x) < 1L || anyNA(value) || any(!is.finite(as.numeric(x))) ||
      any(as.numeric(x) != value) || any(value < minimum)) {
    stop(
      name, " must contain integers greater than or equal to ", minimum, ".",
      call. = FALSE
    )
  }
  value
}


.popart_probability <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || !is.finite(value) || value <= 0 || value >= 1) {
    stop(name, " must be a single number strictly between 0 and 1.", call. = FALSE)
  }
  value
}


.popart_flag <- function(x, name) {
  if (length(x) != 1L || is.na(x) || !is.logical(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}


.validate_popart_control <- function(control) {
  if (!inherits(control, "popart_control")) {
    stop("control must be created by popart_control().", call. = FALSE)
  }
  control
}
