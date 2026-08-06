# Repository-only simulation loader

#' Load the reference-study implementation.
#'
#' The simulation implementation is deliberately kept outside the installed
#' `popart` package. This loader evaluates it in a caller-supplied environment
#' so the research workflow does not add simulation helpers to the package API.
#'
#' @param project_dir Path to the repository root.
#' @param envir Environment that receives the research-only functions.
#'
#' @return Invisibly returns `envir`.
source_reference_study <- function(project_dir = ".", envir = parent.frame()) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  simulation_dir <- file.path(project_dir, "simulation", "R")
  files <- c(
    "configuration.R",
    "persistence.R",
    "data_generation.R",
    "hal_fitting.R",
    "reference_estimators.R",
    "hal_scheduler.R",
    "monte_carlo.R",
    "report_helpers.R",
    "simulation_reporting.R"
  )
  paths <- file.path(simulation_dir, files)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop(
      "Missing simulation source file(s): ",
      paste(basename(missing), collapse = ", "),
      call. = FALSE
    )
  }

  for (path in paths) {
    sys.source(path, envir = envir, chdir = FALSE)
  }
  invisible(envir)
}
