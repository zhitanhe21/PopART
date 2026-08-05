#' Source all functions required by the local AIPW demo.
#'
#' This helper loads the implementation files in dependency order. It is the
#' only source call that user-facing scripts need before running the demo
#' workflow.
#'
#' @param project_dir Character path to the `aipw_local_demo_2` directory.
#'
#' @return Invisibly returns the list produced by `lapply(files, source)`.
source_aipw_local_demo <- function(project_dir) {
  r_dir <- file.path(project_dir, "R")

  files <- c(
    file.path(r_dir, "setup_helpers.R"),
    file.path(r_dir, "persistence.R"),
    file.path(r_dir, "report_helpers.R"),
    file.path(r_dir, "sim2_data.R"),
    file.path(r_dir, "model_fitting.R"),
    file.path(r_dir, "global_scheduler.R"),
    file.path(r_dir, "aipw_estimators.R"),
    file.path(r_dir, "sim2_demo.R"),
    file.path(r_dir, "profiles.R"),
    file.path(r_dir, "sim2_demo_report.R")
  )

  invisible(lapply(files, source))
}

source_fitcvreuse_local <- source_aipw_local_demo
source_aipw_local_demo_2 <- source_aipw_local_demo
