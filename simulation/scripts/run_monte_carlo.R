#!/usr/bin/env Rscript

# Run the repository-only PopART manuscript-reference Monte Carlo study.

`%||%` <- function(x, y) if (is.null(x)) y else x

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

parse_arguments <- function(args) {
  out <- list()
  ii <- 1L
  while (ii <= length(args)) {
    token <- args[[ii]]
    if (!startsWith(token, "--")) {
      stop("Unexpected argument: ", token, call. = FALSE)
    }
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      parts <- strsplit(token, "=", fixed = TRUE)[[1L]]
      out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
      ii <- ii + 1L
      next
    }
    if (token %in% c(
      "help", "dry-run", "no-cache", "no-resume", "no-figures"
    )) {
      out[[token]] <- TRUE
      ii <- ii + 1L
      next
    }
    if (ii == length(args)) {
      stop("Missing value for --", token, call. = FALSE)
    }
    out[[token]] <- args[[ii + 1L]]
    ii <- ii + 2L
  }
  out
}

integer_values <- function(value, name) {
  values <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (length(values) < 1L || anyNA(values) || any(values < 1L)) {
    stop("--", name, " must contain comma-separated positive integers.", call. = FALSE)
  }
  values
}

integer_value <- function(value, name, default = NULL) {
  if (is.null(value)) return(default)
  result <- suppressWarnings(as.integer(value))
  if (length(result) != 1L || is.na(result)) {
    stop("--", name, " must be an integer.", call. = FALSE)
  }
  result
}

print_usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript simulation/scripts/run_monte_carlo.R [options]",
    "",
    "Required sample-size inputs:",
    "  --n-trial 500,5000",
    "  --n-auxiliary 500,5000",
    "",
    "Study options:",
    "  --replicates 20",
    "  --run-id reference_20",
    "  --base-seed 91000000",
    "  --n-clusters 20",
    "  --output-dir simulation/results",
    "  --cache-dir PATH",
    "",
    "Parallel options:",
    "  --n-cores INTEGER",
    "  --n-reserved-cores INTEGER",
    "  --n-cv-workers-per-fit INTEGER",
    "  --max-concurrent-fits INTEGER",
    "  --max-active-replicates INTEGER",
    "  --n-cv-folds INTEGER",
    "  --n-lambda-values INTEGER",
    "",
    "Flags:",
    "  --dry-run --no-cache --no-resume --no-figures",
    sep = "\n"
  ))
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_usage()
  quit(save = "no", status = 0L)
}

project_dir <- find_repository_root()
source(file.path(project_dir, "simulation", "R", "load_simulation.R"))
study <- new.env(parent = globalenv())
source_reference_study(project_dir, envir = study)

n_trial <- integer_values(args[["n-trial"]] %||% "500,5000", "n-trial")
n_auxiliary <- integer_values(
  args[["n-auxiliary"]] %||% "500,5000",
  "n-auxiliary"
)
sample_sizes <- study$sample_size_grid(n_trial, n_auxiliary)

output_dir <- args[["output-dir"]] %||%
  file.path(project_dir, "simulation", "results")
cache_dir <- args[["cache-dir"]] %||%
  file.path(output_dir, "cache")
run_id <- args[["run-id"]] %||%
  paste0("study_", format(Sys.time(), "%Y%m%d_%H%M%S"))

result <- study$run_monte_carlo_study(
  sample_sizes = sample_sizes,
  project_dir = project_dir,
  run_id = run_id,
  base_seed = integer_value(args[["base-seed"]], "base-seed", 91000000L),
  out_dir = output_dir,
  cache_dir = cache_dir,
  use_cache = !isTRUE(args[["no-cache"]]),
  resume = !isTRUE(args[["no-resume"]]),
  m = integer_value(args[["n-clusters"]], "n-clusters", 20L),
  mc_reps = integer_value(args$replicates, "replicates", 1L),
  total_cores = integer_value(args[["n-cores"]], "n-cores"),
  reserve_cores = integer_value(
    args[["n-reserved-cores"]], "n-reserved-cores"
  ),
  cv_workers_per_fit = integer_value(
    args[["n-cv-workers-per-fit"]], "n-cv-workers-per-fit"
  ),
  max_concurrent_hal_fits = integer_value(
    args[["max-concurrent-fits"]], "max-concurrent-fits"
  ),
  max_active_mc = integer_value(
    args[["max-active-replicates"]], "max-active-replicates"
  ),
  cv_folds = integer_value(args[["n-cv-folds"]], "n-cv-folds"),
  cv_nlambda = integer_value(
    args[["n-lambda-values"]], "n-lambda-values"
  ),
  dry_run = isTRUE(args[["dry-run"]])
)

if (!isTRUE(args[["dry-run"]])) {
  study$write_monte_carlo_report(run_id = run_id, out_dir = output_dir)
  if (!isTRUE(args[["no-figures"]])) {
    study$write_monte_carlo_figures(run_id = run_id, out_dir = output_dir)
  }
}

invisible(result)
