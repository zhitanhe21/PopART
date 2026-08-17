###############################################################################
# Reference-study configuration and parallel planning helpers
###############################################################################


# Shared statistical defaults for the experimental laptop configuration.
# Keep these values centralized because they also participate in cache and
# checkpoint identities.
REFERENCE_DEFAULT_CV_FOLDS <- 3L
REFERENCE_DEFAULT_NUM_KNOTS <- c(5, 3)


#' Load packages required by the reference study.
#'
#' This helper centralizes package loading for the research scripts and functions.
#' Keeping package loading in one place makes the direct scripts shorter and
#' makes missing-package errors easier to diagnose.
#'
#' @return Invisibly returns `NULL`. In package mode, required namespaces are
#'   loaded without changing the search path. The legacy source mode also
#'   attaches them so unqualified source-file calls keep working.
load_reference_study_packages <- function() {
  packages <- c(
    "dplyr", "tidyr", "ggplot2", "hal9001", "doParallel", "parallel"
  )
  namespaces <- lapply(packages, loadNamespace)

  # Functions sourced into the global environment do not receive package
  # namespace imports. Preserve that legacy entry point without attaching
  # dependencies when this helper is called from the installed package.
  package_mode <- isNamespace(environment(load_reference_study_packages))
  if (!isTRUE(package_mode)) {
    for (i in seq_along(packages)) {
      search_name <- paste0("package:", packages[[i]])
      if (!search_name %in% search()) {
        suppressPackageStartupMessages(attachNamespace(namespaces[[i]]))
      }
    }
  }

  invisible(NULL)
}


#' Limit threaded math libraries to one thread.
#'
#' The reference study uses explicit R-level parallelism for HAL fits and cross-validation.
#' This helper prevents BLAS/OpenMP-style libraries from silently creating
#' additional nested threads.
#'
#' @return Invisibly returns the result of `Sys.setenv()`.
set_single_threaded_math <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
}


#' Parse command-line flags into a named list.
#'
#' Supports flags in either `--key value` or `--key=value` form. Dashes in flag
#' names are converted to underscores, so `--cv-nlambda 5` becomes
#' `args$cv_nlambda`.
#'
#' @param args Character vector, usually `commandArgs(trailingOnly = TRUE)`.
#'
#' @return A named list of parsed argument values. Standalone flags are stored as
#'   `TRUE`.
parse_cli_args <- function(args) {
  out <- list()
  ii <- 1L
  while (ii <= length(args)) {
    arg <- args[[ii]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }

    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", key_value[[1]])
    if (length(key_value) > 1L) {
      out[[key]] <- paste(key_value[-1L], collapse = "=")
      ii <- ii + 1L
    } else if (ii == length(args) || startsWith(args[[ii + 1L]], "--")) {
      out[[key]] <- TRUE
      ii <- ii + 1L
    } else {
      out[[key]] <- args[[ii + 1L]]
      ii <- ii + 2L
    }
  }
  out
}


arg_value <- function(args, name, default = NULL) {
  if (!is.null(args[[name]])) args[[name]] else default
}


arg_int <- function(args, name, default = NA_integer_) {
  value <- arg_value(args, name, default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) default else out
}


arg_bool <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, default)
  if (is.logical(value)) return(value)
  tolower(as.character(value)) %in% c("1", "true", "t", "yes", "y")
}


get_project_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
  } else {
    script_dir <- getwd()
  }
  normalizePath(file.path(script_dir, ".."))
}


set_int_env_from_arg <- function(args, arg_name, env_name) {
  value <- arg_int(args, arg_name, NA_integer_)
  if (!is.na(value)) {
    do.call(Sys.setenv, as.list(stats::setNames(as.character(value), env_name)))
  }
  invisible(value)
}


read_positive_int_env <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) default else out
}


read_nonnegative_int_env <- function(name, default = NA_integer_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 0L) default else out
}


read_positive_numeric_env <- function(name, default = NA_real_) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out) || out <= 0) default else out
}


#' Detect the local CPU budget.
#'
#' CPU detection respects, in order, `HAL_TOTAL_CORES`, `SLURM_CPUS_PER_TASK`,
#' and local core detection. The returned value is used only for planning; it
#' does not allocate cores by itself.
#'
#' @param prefer_physical Logical, whether to prefer physical cores when R can
#'   detect both physical and logical cores. The default is `FALSE` because the
#'   local worker budget is expressed in available logical CPU slots.
#'
#' @return A positive integer CPU count.
#' @export
detect_local_cores <- function(prefer_physical = FALSE) {
  override <- read_positive_int_env("HAL_TOTAL_CORES", NA_integer_)
  if (!is.na(override)) return(override)

  slurm_cores <- read_positive_int_env("SLURM_CPUS_PER_TASK", NA_integer_)
  if (!is.na(slurm_cores)) return(slurm_cores)

  available <- tryCatch(
    suppressWarnings(as.integer(parallelly::availableCores()[[1L]])),
    error = function(e) NA_integer_
  )
  nproc <- tryCatch(
    suppressWarnings(as.integer(system2("nproc", stdout = TRUE, stderr = FALSE)[[1]])),
    error = function(e) NA_integer_
  )

  physical <- suppressWarnings(parallel::detectCores(logical = FALSE))
  logical <- suppressWarnings(parallel::detectCores(logical = TRUE))
  detected <- if (prefer_physical && !is.na(physical)) {
    physical
  } else if (!is.na(available)) {
    available
  } else {
    logical
  }
  if (is.na(detected)) detected <- 1L

  if (!is.na(available)) detected <- min(detected, available)
  if (!is.na(nproc)) detected <- min(detected, nproc)
  max(1L, as.integer(detected))
}


#' Choose a local global-HAL scheduling plan.
#'
#' The default uses serial CV inside each complete HAL fit, reserves two
#' available logical CPU slots, lets each run occupy at most two HAL slots, and
#' admits 20 percent more run contexts than the minimum needed to feed all
#' global slots. The active-run value is a memory valve, not another
#' multiplicative CPU layer.
#'
#' @param total_cores Integer local CPU budget.
#' @param reserve_cores Integer number of cores to leave unused.
#' @param n_fit_jobs Integer number of nuisance-model fit jobs.
#' @param cv_folds Integer number of HAL CV folds.
#' @param n_runs Integer number of sample-combination-by-MC runs.
#' @param active_mc_headroom Positive scheduling headroom multiplier.
#'
#' @return A list with `total_cores`, `reserve_cores`, `usable_cores`,
#'   per-run and global HAL limits, and the active-run window.
#' @export
choose_fitcv_plan <- function(total_cores = detect_local_cores(),
                              reserve_cores = 2L,
                              n_fit_jobs = 4L,
                              cv_folds = REFERENCE_DEFAULT_CV_FOLDS,
                              n_runs = 1L,
                              active_mc_headroom = 1.2) {
  total_cores <- suppressWarnings(as.integer(total_cores))
  if (is.na(total_cores) || total_cores < 1L) total_cores <- detect_local_cores()
  reserve_cores <- suppressWarnings(as.integer(reserve_cores))
  if (is.na(reserve_cores) || reserve_cores < 0L) reserve_cores <- 2L
  n_fit_jobs <- suppressWarnings(as.integer(n_fit_jobs))
  if (is.na(n_fit_jobs) || n_fit_jobs < 1L) n_fit_jobs <- 1L
  n_runs <- suppressWarnings(as.integer(n_runs))
  if (is.na(n_runs) || n_runs < 1L) n_runs <- 1L
  cv_folds <- suppressWarnings(as.integer(cv_folds))
  if (is.na(cv_folds) || cv_folds < 1L) cv_folds <- 1L
  active_mc_headroom <- suppressWarnings(as.numeric(active_mc_headroom))
  if (!is.finite(active_mc_headroom) || active_mc_headroom <= 0) {
    active_mc_headroom <- 1.2
  }
  reserve_cores <- if (total_cores <= 1L) {
    0L
  } else {
    min(reserve_cores, total_cores - 1L)
  }
  usable_cores <- max(1L, total_cores - reserve_cores)

  cv_workers_per_fit <- 1L
  fit_workers <- min(2L, usable_cores, n_fit_jobs)
  max_concurrent_hal_fits <- min(
    max(1L, floor(usable_cores / cv_workers_per_fit)),
    n_runs * fit_workers
  )
  max_active_mc <- min(
    n_runs,
    max(
      1L,
      as.integer(ceiling(
        active_mc_headroom * max_concurrent_hal_fits / fit_workers
      ))
    )
  )

  list(
    total_cores = total_cores,
    reserve_cores = reserve_cores,
    usable_cores = usable_cores,
    fit_workers = as.integer(fit_workers),
    fit_workers_per_mc = as.integer(fit_workers),
    cv_workers_per_fit = as.integer(cv_workers_per_fit),
    max_concurrent_hal_fits = as.integer(max_concurrent_hal_fits),
    global_hal_slots = as.integer(max_concurrent_hal_fits),
    max_active_mc = as.integer(max_active_mc),
    active_mc_headroom = active_mc_headroom
  )
}


#' Build the local parallel configuration.
#'
#' Reads environment overrides such as `HAL_TOTAL_CORES`, `HAL_FIT_WORKERS`,
#' `HAL_MAX_CONCURRENT_HAL_FITS`, `HAL_MAX_ACTIVE_MC`, and the HAL CV controls,
#' then returns the configuration used by the global scheduler. Manual values
#' are clipped so the complete-fit slots stay within the usable CPU budget.
#'
#' @param n_fit_jobs Integer number of nuisance-model fit jobs.
#' @param default_cv_folds Integer default HAL CV folds.
#' @param default_cv_nlambda Integer default number of HAL lambda values.
#' @param backend Character dispatch backend for fit-level parallelism:
#'   `"auto"`, `"fork"`, or `"psock"`.
#' @param n_runs Integer number of sample-combination-by-MC runs.
#' @param default_active_mc_headroom Positive automatic active-window
#'   multiplier. The default `1.2` matches the adopted global-queue rule.
#'
#' @return A named list with CPU, worker, CV, and backend settings.
#' @export
configure_fitcv_parallel <- function(n_fit_jobs = 4L,
                                     default_cv_folds = REFERENCE_DEFAULT_CV_FOLDS,
                                     default_cv_nlambda = 50L,
                                     backend = Sys.getenv("HAL_FIT_BACKEND", "auto"),
                                     n_runs = 1L,
                                     default_active_mc_headroom = 1.2) {
  total_cores <- detect_local_cores()
  cv_folds <- read_positive_int_env("HAL_CV_FOLDS", default_cv_folds)
  cv_nlambda <- read_positive_int_env("HAL_CV_NLAMBDA", default_cv_nlambda)
  if (cv_folds < 3L) {
    stop("HAL_CV_FOLDS/cv_folds must be at least 3.", call. = FALSE)
  }
  reserve_cores <- read_nonnegative_int_env("HAL_RESERVE_CORES", 2L)
  active_mc_headroom <- read_positive_numeric_env(
    "HAL_ACTIVE_MC_HEADROOM",
    default_active_mc_headroom
  )
  plan <- choose_fitcv_plan(
    total_cores = total_cores,
    reserve_cores = reserve_cores,
    n_fit_jobs = n_fit_jobs,
    cv_folds = cv_folds,
    n_runs = n_runs,
    active_mc_headroom = active_mc_headroom
  )

  cv_workers_per_fit <- read_positive_int_env(
    "HAL_CV_WORKERS_PER_FIT",
    plan$cv_workers_per_fit
  )
  cv_workers_per_fit <- min(max(1L, cv_workers_per_fit), cv_folds)
  cv_workers_per_fit <- min(cv_workers_per_fit, plan$usable_cores)

  cpu_hal_slots <- as.integer(max(
    1L,
    floor(plan$usable_cores / cv_workers_per_fit)
  ))
  fit_workers <- read_positive_int_env("HAL_FIT_WORKERS", 2L)
  fit_workers <- min(max(1L, fit_workers), n_fit_jobs, cpu_hal_slots)

  max_concurrent_hal_fits <- read_positive_int_env(
    "HAL_MAX_CONCURRENT_HAL_FITS",
    cpu_hal_slots
  )
  max_concurrent_hal_fits <- as.integer(min(
    max(1L, max_concurrent_hal_fits),
    cpu_hal_slots,
    max(1L, as.integer(n_runs)) * fit_workers
  ))

  auto_active_mc <- min(
    max(1L, as.integer(n_runs)),
    max(
      1L,
      as.integer(ceiling(
        active_mc_headroom * max_concurrent_hal_fits / fit_workers
      ))
    )
  )
  max_active_mc <- read_positive_int_env("HAL_MAX_ACTIVE_MC", auto_active_mc)
  max_active_mc <- as.integer(min(
    max(1L, max_active_mc),
    max(1L, as.integer(n_runs))
  ))
  max_concurrent_hal_fits <- as.integer(min(
    max_concurrent_hal_fits,
    max_active_mc * fit_workers
  ))

  list(
    total_cores = total_cores,
    reserve_cores = plan$reserve_cores,
    usable_cores = plan$usable_cores,
    fit_workers = fit_workers,
    fit_workers_per_mc = fit_workers,
    cv_workers_per_fit = cv_workers_per_fit,
    max_concurrent_hal_fits = as.integer(max_concurrent_hal_fits),
    global_hal_slots = as.integer(max_concurrent_hal_fits),
    max_active_mc = as.integer(max_active_mc),
    active_mc_headroom = active_mc_headroom,
    cv_folds = cv_folds,
    cv_nlambda = cv_nlambda,
    backend = backend
  )
}


#' Create an explicit sample-size grid for the reference study.
#'
#' @param n_trial Integer vector of trial sample sizes.
#' @param n_auxiliary Integer vector of auxiliary sample sizes.
#'
#' @return A data frame with columns `n_trial` and `n_auxiliary`.
sample_size_grid <- function(n_trial, n_auxiliary) {
  n_trial <- as.integer(n_trial)
  n_auxiliary <- as.integer(n_auxiliary)
  if (length(n_trial) < 1L || length(n_auxiliary) < 1L ||
      anyNA(n_trial) || anyNA(n_auxiliary) ||
      any(n_trial < 1L) || any(n_auxiliary < 1L)) {
    stop("n_trial and n_auxiliary must contain positive integers.", call. = FALSE)
  }
  expand.grid(
    n_trial = n_trial,
    n_auxiliary = n_auxiliary,
    KEEP.OUT.ATTRS = FALSE
  )
}


#' Expand explicit sample-size settings into a Monte Carlo grid.
#'
#' Each Monte Carlo replicate uses one seed across all requested sample-size
#' combinations, matching the manuscript reference-study seed organization.
#'
#' @param sizes Data frame with `n_trial` and `n_aux` columns.
#' @param mc_reps Positive integer number of Monte Carlo replicates per
#'   sample-size combination.
#' @param base_seed Nonnegative integer seed offset.
#'
#' @return A data frame with the original size columns plus `mc_rep`,
#'   `mc_reps`, and `seed`.
#' @export
expand_monte_carlo_grid <- function(sizes, mc_reps = 1L,
                                    base_seed = 91000000L) {
  if (!is.data.frame(sizes) ||
      !all(c("n_trial", "n_aux") %in% names(sizes)) ||
      nrow(sizes) < 1L) {
    stop(
      "sizes must be a nonempty data frame with n_trial and n_aux.",
      call. = FALSE
    )
  }

  mc_reps <- suppressWarnings(as.integer(mc_reps))
  base_seed <- suppressWarnings(as.integer(base_seed))
  if (is.na(mc_reps) || mc_reps < 1L) {
    stop("mc_reps must be a positive integer.", call. = FALSE)
  }
  if (is.na(base_seed) || base_seed < 0L ||
      base_seed > .Machine$integer.max - mc_reps) {
    stop(
      "base_seed must leave room for all Monte Carlo replicate seeds.",
      call. = FALSE
    )
  }

  n_sizes <- nrow(sizes)
  out <- sizes[
    rep(seq_len(n_sizes), times = mc_reps),
    ,
    drop = FALSE
  ]
  out$mc_rep <- rep(seq_len(mc_reps), each = n_sizes)
  out$mc_reps <- mc_reps
  out$seed <- base_seed + out$mc_rep
  rownames(out) <- NULL
  out
}


#' Print the selected local parallel plan.
#'
#' @param config List returned by `configure_fitcv_parallel()`.
#'
#' @return Invisibly returns `NULL`; called for console output.
print_parallel_plan <- function(config) {
  cat("Global HAL scheduler plan\n")
  cat("  total_cores = ", config$total_cores, "\n", sep = "")
  cat("  reserve_cores = ", config$reserve_cores, "\n", sep = "")
  cat("  usable_cores = ", config$usable_cores, "\n", sep = "")
  cat("  fit_workers_per_mc = ", config$fit_workers, "\n", sep = "")
  cat("  cv_workers_per_fit = ", config$cv_workers_per_fit, "\n", sep = "")
  cat("  max_concurrent_hal_fits = ",
      config$max_concurrent_hal_fits, "\n", sep = "")
  cat("  max_active_mc = ", config$max_active_mc, "\n", sep = "")
  cat("  active_mc_headroom = ",
      config$active_mc_headroom, "\n", sep = "")
  cat("  max_requested_cpu_units = ",
      config$max_concurrent_hal_fits * config$cv_workers_per_fit,
      "\n", sep = "")
  cat("  cv_folds = ", config$cv_folds, "\n", sep = "")
  cat("  cv_nlambda = ", config$cv_nlambda, "\n", sep = "")
  if (!is.null(config$scheduler_backend)) {
    cat("  scheduler_backend = ", config$scheduler_backend, "\n", sep = "")
  }
  cat("  legacy_fit_backend = ", config$backend, "\n", sep = "")
}
