#!/usr/bin/env Rscript

# Regenerate a report and figures from a saved Monte Carlo study.

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
    key <- args[[ii]]
    if (!startsWith(key, "--") || ii == length(args)) {
      stop("Arguments must use --name value pairs.", call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[ii + 1L]]
    ii <- ii + 2L
  }
  out
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))
run_id <- args[["run-id"]]
if (is.null(run_id) || !nzchar(run_id)) {
  stop("--run-id is required.", call. = FALSE)
}

project_dir <- find_repository_root()
source(file.path(project_dir, "simulation", "R", "load_simulation.R"))
study <- new.env(parent = globalenv())
source_reference_study(project_dir, envir = study)
study$load_reference_study_packages()

output_dir <- args[["output-dir"]] %||%
  file.path(project_dir, "simulation", "results")
figure_dir <- args[["figure-dir"]] %||% file.path(output_dir, "figures")

study$write_monte_carlo_report(run_id = run_id, out_dir = output_dir)
study$write_monte_carlo_figures(
  run_id = run_id,
  out_dir = output_dir,
  figure_dir = figure_dir
)
