###############################################################################
# User-facing estimation controls
###############################################################################


#' Configure PopART estimation
#'
#' Creates and validates the computational and inferential settings used by
#' [fit_popart()]. The settings control the four highly adaptive lasso (HAL)
#' nuisance fits, reproducible cross-validation folds, Wald intervals, weight
#' normalization, and diagnostic reporting.
#'
#' Parallelism is available at one level at a time. Set `n_fit_workers > 1` to
#' run complete nuisance fits concurrently, or set `n_cv_workers > 1` to run
#' folds concurrently while the four nuisance fits are serial. A worker is an R
#' process scheduled by the operating system; it is not a promise of a dedicated
#' physical CPU core.
#'
#' By default, `n_fit_workers = "auto"` queries the CPU slots available to the R
#' process with [parallelly::availableCores()]. The resolved value is
#' `max(1, min(4, available_cpu_slots - 2))`. This leaves at least two detected
#' CPU slots outside the complete-fit pool when possible and caps the result at four
#' because one PopART analysis has four nuisance fits.
#' If `n_cv_workers > 1`, automatic fit-level parallelism resolves to one so
#' that only the CV folds run concurrently.
#'
#' `parallel_backend = "auto"` selects PSOCK on Windows and forked processes on
#' other supported operating systems. The HAL basis arguments are passed to
#' `hal9001::fit_hal()`. Auxiliary-weight normalization rescales the supplied
#' weights relative to unit-weighted trial records and records the scale factor
#' in the fitted object's diagnostics.
#'
#' @param n_cv_folds A single integer greater than or equal to 3 giving the
#'   number of cross-validation folds for each HAL fit. Default: `3L`.
#' @param n_lambda_values A single positive integer giving the number of lasso
#'   penalty values evaluated within each HAL fit. Default: `50L`.
#' @param n_cv_workers A single positive integer, no greater than `n_cv_folds`,
#'   giving the number of R workers used across cross-validation folds. It must
#'   equal one when a manually specified `n_fit_workers` is greater than one.
#'   Default: `1L`.
#' @param n_fit_workers Either `"auto"` or a single integer from 1 through 4
#'   giving the maximum number of complete nuisance fits run concurrently.
#'   Automatic mode uses the conservative CPU-based rule described above and
#'   resolves to one when `n_cv_workers > 1`. A manually specified value greater
#'   than one cannot be combined with `n_cv_workers > 1`. Default: `"auto"`.
#' @param parallel_backend A character string selecting the complete-fit
#'   parallel backend: `"auto"`, `"fork"`, or `"psock"`. The fork backend is
#'   unavailable on Windows. Default: `"auto"`.
#' @param random_seed A single nonnegative integer used to construct reproducible
#'   cross-validation fold assignments. Default: `1L`.
#' @param smoothness_orders A nonempty vector of nonnegative integers specifying
#'   HAL basis smoothness orders. Default: `1L`.
#' @param max_degree A single positive integer giving the maximum HAL interaction
#'   degree. Default: `2L`.
#' @param num_knots A nonempty vector of positive integers specifying the HAL
#'   knot counts. Default: `c(5L, 3L)`.
#' @param treatment_probability A single numeric value strictly between zero and
#'   one giving the known probability of assignment to the treated arm. Default:
#'   `0.5`, corresponding to equal randomization.
#' @param conf_level A single numeric value strictly between zero and one giving
#'   the confidence level for Wald intervals. Default: `0.95`.
#' @param normalize_auxiliary_weights A single logical value indicating whether
#'   auxiliary weights are divided by their mean before fitting. Default: `TRUE`.
#' @param keep_nuisance_fits A single logical value indicating whether the four
#'   fitted HAL objects, including their underlying `cv.glmnet` fits, are
#'   retained in the returned `popart_fit` object. When `FALSE`, HAL still runs
#'   the same cross-validation and produces the same predictions, but its
#'   unneeded `cv.glmnet` object is discarded to reduce memory use and worker
#'   transfer overhead. This does not skip `glmnet` computation or reduce the
#'   peak memory used while an individual model is being fitted. Default:
#'   `FALSE`.
#' @param positivity_threshold A single numeric value strictly between zero and
#'   `0.5`. Estimated probabilities closer than this value to a boundary are
#'   flagged in diagnostics but are not truncated. Default: `0.01`.
#'
#' @return An object of class `popart_control`, implemented as a list with the
#'   following elements:
#'   \itemize{
#'   \item `n_cv_folds`, `n_lambda_values`, `n_cv_workers`, and `n_fit_workers`:
#'     validated integer fitting and resolved worker counts.
#'   \item `n_fit_workers_auto`: whether the complete-fit worker count was
#'     selected automatically.
#'   \item `detected_cpu_slots`: the positive number of logical CPU slots or
#'     scheduler-limited slots reported as available when the control object was
#'     created.
#'   \item `parallel_backend`: the requested backend string; `"auto"` is
#'     resolved when fitting begins.
#'   \item `random_seed`: the validated integer cross-validation seed.
#'   \item `smoothness_orders`, `max_degree`, and `num_knots`: validated HAL
#'     basis settings.
#'   \item `treatment_probability`: the treated-arm randomization probability.
#'   \item `conf_level`: the Wald confidence level.
#'   \item `normalize_auxiliary_weights`: the weight-normalization flag.
#'   \item `keep_nuisance_fits`: the nuisance-fit retention flag.
#'   \item `positivity_threshold`: the diagnostic probability threshold.
#'   }
#'
#' @examples
#' # The default selects a conservative complete-fit worker count automatically.
#' automatic <- popart_control()
#' automatic
#'
#' # Request fully serial execution explicitly.
#' serial <- popart_control(n_fit_workers = 1L, n_cv_workers = 1L)
#'
#' # Use three workers for cross-validation folds.
#' cv_parallel <- popart_control(n_fit_workers = 1L, n_cv_workers = 3L)
#'
#' # Or run two of the four complete nuisance fits concurrently.
#' fit_parallel <- popart_control(n_fit_workers = 2L, n_cv_workers = 1L)
#'
#' @export
popart_control <- function(
    n_cv_folds = 3L,
    n_lambda_values = 50L,
    n_cv_workers = 1L,
    n_fit_workers = "auto",
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
  detected_cpu_slots <- .popart_available_cpu_slots()
  n_fit_workers_auto <- is.character(n_fit_workers) &&
    length(n_fit_workers) == 1L &&
    !is.na(n_fit_workers) &&
    identical(n_fit_workers, "auto")
  if (n_fit_workers_auto) {
    n_fit_workers <- if (n_cv_workers > 1L) {
      1L
    } else {
      .popart_auto_fit_workers(detected_cpu_slots)
    }
  } else {
    if (is.character(n_fit_workers)) {
      stop(
        "n_fit_workers must be \"auto\" or an integer between 1 and 4.",
        call. = FALSE
      )
    }
    n_fit_workers <- .popart_count(
      n_fit_workers,
      "n_fit_workers",
      maximum = 4L
    )
  }
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
      n_fit_workers_auto = n_fit_workers_auto,
      detected_cpu_slots = detected_cpu_slots,
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


#' Print PopART estimation controls
#'
#' Displays the principal cross-validation, parallelism, HAL-basis, and
#' confidence-level settings in a `popart_control` object.
#'
#' @param x A `popart_control` object created by [popart_control()].
#' @param ... Additional arguments, currently ignored.
#'
#' @return The input `popart_control` object, returned invisibly.
#'
#' @examples
#' print(popart_control())
#'
#' @export
print.popart_control <- function(x, ...) {
  cat("PopART estimation controls\n")
  cat("  CV folds / workers: ", x$n_cv_folds, " / ", x$n_cv_workers, "\n", sep = "")
  worker_note <- if (isTRUE(x$n_fit_workers_auto)) {
    if (x$n_cv_workers > 1L) {
      paste0(
        " (auto: fit-level parallelism disabled while CV uses ",
        x$n_cv_workers,
        " workers)"
      )
    } else {
      paste0(" (auto from ", x$detected_cpu_slots, " available CPU slots)")
    }
  } else {
    " (manual)"
  }
  cat(
    "  Complete-fit workers: ",
    x$n_fit_workers,
    worker_note,
    "\n",
    sep = ""
  )
  cat("  Parallel backend: ", x$parallel_backend, "\n", sep = "")
  cat("  HAL knots: ", paste(x$num_knots, collapse = ", "), "\n", sep = "")
  cat("  Confidence level: ", x$conf_level, "\n", sep = "")
  invisible(x)
}


.popart_available_cpu_slots <- function() {
  available <- tryCatch(
    parallelly::availableCores(),
    error = function(condition) 1L
  )
  if (length(available) < 1L) {
    return(1L)
  }
  available <- suppressWarnings(as.integer(available[[1L]]))
  if (length(available) != 1L || is.na(available) || available < 1L) {
    return(1L)
  }
  available
}


.popart_auto_fit_workers <- function(available_cpu_slots) {
  numeric_cpu_slots <- suppressWarnings(as.numeric(available_cpu_slots))
  available_cpu_slots <- suppressWarnings(as.integer(available_cpu_slots))
  if (length(numeric_cpu_slots) != 1L || !is.finite(numeric_cpu_slots) ||
      length(available_cpu_slots) != 1L || is.na(available_cpu_slots) ||
      numeric_cpu_slots != available_cpu_slots || available_cpu_slots < 1L) {
    stop("available_cpu_slots must be a positive integer.", call. = FALSE)
  }
  min(4L, max(1L, available_cpu_slots - 2L))
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
