###############################################################################
###############################################################################

# Local optimized PopART Simulation 2 demo v2 analysis

###############################################################################
###############################################################################


# setup -------------------------------------------------------------------

rm(list = ls())

find_script_dir <- function() {
  args_all <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_all, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(path)) {
      return(dirname(normalizePath(path)))
    }
  }

  getwd()
}

find_project_dir <- function(script_dir) {
  raw_candidates <- c(
    script_dir,
    file.path(script_dir, ".."),
    file.path(script_dir, "..", ".."),
    getwd(),
    file.path(getwd(), "aipw_local_demo_2")
  )

  candidates <- unique(normalizePath(raw_candidates, mustWork = FALSE))
  for (path in candidates) {
    current <- path
    for (ii in seq_len(5L)) {
      if (file.exists(file.path(current, "R", "source_all.R")) &&
          file.exists(file.path(current, "simulation", "sim_analysis", "sim2_demo_analysis.R"))) {
        return(normalizePath(current))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  stop(
    "Could not locate aipw_local_demo_2. Open this script from inside the demo ",
    "folder, or set the RStudio working directory to the HAL project root.",
    call. = FALSE
  )
}

launch_dir <- normalizePath(getwd())
script_dir <- find_script_dir()
project_dir <- find_project_dir(script_dir)
setwd(project_dir)

source(file.path(project_dir, "R", "source_all.R"))
source_aipw_local_demo(project_dir)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

resolve_cli_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(path)
  is_absolute <- grepl(
    "^(?:[A-Za-z]:[/\\\\]|[/\\\\]{2}|/)",
    path,
    perl = TRUE
  )
  candidate <- if (is_absolute) path else file.path(launch_dir, path)
  normalizePath(candidate, winslash = "/", mustWork = FALSE)
}

allowed_args <- c(
  "help", "cluster_id", "scale", "run_id", "out_dir", "figure_dir",
  "no_plots"
)
unknown_args <- setdiff(names(args), allowed_args)
if (length(unknown_args) > 0L) {
  stop(
    "Unknown option(s): --",
    paste(gsub("_", "-", unknown_args), collapse = ", --"),
    ". Use --help to see supported options.",
    call. = FALSE
  )
}
if (arg_bool(args, "help", FALSE)) {
  cat(paste(
    "AIPW Local Demo 2 analysis",
    "",
    "Usage:",
    "  Rscript simulation/sim_analysis/sim2_demo_analysis.R [options]",
    "",
    "Options:",
    "  --cluster-id N",
    "  --scale demo|small|large",
    "  --run-id ID",
    "  --out-dir PATH",
    "  --figure-dir PATH",
    "  --no-plots",
    "  --help",
    sep = "\n"
  ))
  cat("\n")
  quit(save = "no", status = 0L)
}


# analysis parameters -----------------------------------------------------

cluster_id <- arg_int(args, "cluster_id", 0L)
scale <- arg_value(args, "scale", "demo")
run_id <- arg_value(args, "run_id", paste0("sd", cluster_id))
out_dir <- resolve_cli_path(arg_value(
  args,
  "out_dir",
  file.path(project_dir, "simulation", "sim_data", "sim2")
))
figure_dir <- resolve_cli_path(arg_value(
  args,
  "figure_dir",
  file.path(project_dir, "simulation", "sim_figures", "sim2")
))
make_plots <- !arg_bool(args, "no_plots", FALSE)

if (is.na(cluster_id) || cluster_id < 0L) {
  stop("--cluster-id must be a nonnegative integer.", call. = FALSE)
}
if (!scale %in% c("demo", "small", "large")) {
  stop("--scale must be one of: demo, small, large.", call. = FALSE)
}
if (!grepl("^[A-Za-z0-9_.-]+$", run_id)) {
  stop(
    "--run-id may contain only letters, digits, underscore, dot, and hyphen.",
    call. = FALSE
  )
}


# summarize simulations ---------------------------------------------------

write_sim2_demo_report(
  project_dir = project_dir,
  scale = scale,
  run_id = run_id,
  out_dir = out_dir
)

if (isTRUE(make_plots)) {
  write_sim2_demo_figures(
    project_dir = project_dir,
    scale = scale,
    run_id = run_id,
    out_dir = out_dir,
    figure_dir = figure_dir
  )
}
