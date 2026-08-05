###############################################################################
###############################################################################

# Local optimized PopART Simulation 2 demo v2 script

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
          file.exists(file.path(current, "simulation", "sim_scripts", "ss2_demo.R"))) {
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

print_cli_help <- function() {
  cat(paste(
    "AIPW Local Demo 2",
    "",
    "Usage:",
    "  Rscript simulation/sim_scripts/ss2_demo.R [options]",
    "",
    "Workload options:",
    "  --profile smoke|quick|hour|formal",
    "  --scale demo|small|large",
    "  --max-runs N  (sample-size combinations before MC expansion)",
    "  --mc-reps N   (replicates per selected sample-size combination)",
    "  --n-trial N --n-aux N",
    "  --run-id ID --base-seed N",
    "  --p-resp 0.5 --p-cens 0.3  (fixed by the original DGP)",
    "  --cv-folds N --cv-nlambda N",
    "",
    "Global HAL scheduler:",
    "  R = selected sample-size-combination x MC runs",
    "  U = usable CPU units after reservation",
    "  c = --cv-workers-per-fit N                 (default 1)",
    "  F = --fit-workers N, the per-run fit cap  (default 2)",
    "  S = floor(U / c), the global HAL-slot CPU ceiling",
    "  A = min(R, ceiling(1.2 * S / F)), the prepared-run window",
    "  --total-cores N --reserve-cores N",
    "     no-profile reserve default: 1; named profiles: 0",
    "  --max-concurrent-hal-fits N  (cap effective S; cannot exceed CPU limit)",
    "  --max-active-mc N            (memory valve; lower uses less RAM)",
    "  --active-mc-headroom X       (replace 1.2 in automatic A)",
    "  --backend auto|psock|fork    (compatibility field; slots use callr)",
    "",
    "  Persistent callr slots share one eligible-fit queue across sample",
    "  sizes and MC replicates. After a run's four fits succeed, AIPW is",
    "  computed and checkpointed immediately and the run context is freed.",
    "  With c > 1, every slot owns a private reusable CV PSOCK pool.",
    "  The callr package is required.",
    "",
    "Persistence/output options:",
    "  --use-cache | --no-cache",
    "  --no-resume | --checkpoint-dir PATH",
    "  --no-analysis | --no-plots",
    "  --dry-run",
    "  --help",
    "",
    "Timing and memory:",
    "  Use invocation_wall_elapsed_sec for end-to-end scheduler runtime.",
    "  Summed per-run elapsed_sec values overlap and are not global wall time.",
    "  Active contexts retain data, four fit inputs, and completed HAL fits;",
    "  lower --max-active-mc, then --max-concurrent-hal-fits, if RAM is tight.",
    sep = "\n"
  ))
  cat("\n")
}

allowed_args <- c(
  "help", "profile", "cluster_id", "scale", "run_id", "base_seed",
  "m", "p_resp", "p_cens", "max_runs", "mc_reps", "n_trial", "n_aux",
  "total_cores", "reserve_cores", "fit_workers", "cv_workers_per_fit",
  "max_concurrent_hal_fits", "max_active_mc", "active_mc_headroom",
  "cv_folds", "cv_nlambda", "out_dir", "cache_dir", "backend",
  "use_cache", "no_cache", "no_resume", "checkpoint_dir", "dry_run",
  "no_analysis", "figure_dir", "no_plots"
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
  print_cli_help()
  quit(save = "no", status = 0L)
}

profile_name <- arg_value(args, "profile", NULL)
profile_defaults <- if (is.null(profile_name)) {
  list()
} else {
  profile_name <- normalize_sim2_profile(profile_name)
  sim2_profile_defaults(profile_name)
}

validate_integer_arg <- function(name, minimum = NULL) {
  value <- args[[name]]
  if (is.null(value)) return(invisible(NULL))
  text <- as.character(value)
  if (length(text) != 1L || !grepl("^-?[0-9]+$", text)) {
    stop("--", gsub("_", "-", name), " must be an integer.", call. = FALSE)
  }
  number <- as.integer(text)
  if (is.na(number) || (!is.null(minimum) && number < minimum)) {
    stop(
      "--", gsub("_", "-", name), " must be >= ", minimum, ".",
      call. = FALSE
    )
  }
  invisible(number)
}

for (name in c("cluster_id", "reserve_cores")) {
  validate_integer_arg(name, 0L)
}
validate_integer_arg("max_runs", 1L)
for (name in c(
  "m", "mc_reps", "n_trial", "n_aux", "total_cores", "fit_workers",
  "cv_workers_per_fit", "max_concurrent_hal_fits", "max_active_mc",
  "cv_nlambda"
)) {
  validate_integer_arg(name, 1L)
}
validate_integer_arg("cv_folds", 3L)
validate_integer_arg("base_seed", 0L)
if (!is.null(args$active_mc_headroom)) {
  active_mc_headroom_check <- suppressWarnings(as.numeric(
    args$active_mc_headroom
  ))
  if (length(active_mc_headroom_check) != 1L ||
      !is.finite(active_mc_headroom_check) ||
      active_mc_headroom_check <= 0) {
    stop("--active-mc-headroom must be a positive number.", call. = FALSE)
  }
}


# simulation parameters ---------------------------------------------------

# These defaults are intentionally small so the script can be run locally.
# The automatic global plan reserves one CPU unit without a profile, uses a
# per-run fit cap of two, and keeps CV serial inside each fit by default.
cluster_id <- arg_int(args, "cluster_id", 0L)
scale <- arg_value(
  args,
  "scale",
  sim2_profile_default(profile_defaults, "scale", "demo")
)
default_run_id <- if (is.null(profile_name)) {
  paste0("sd", cluster_id)
} else {
  paste0("sd", cluster_id, "_", profile_name)
}
run_id <- arg_value(args, "run_id", default_run_id)
base_seed <- arg_int(args, "base_seed", 91000000L + 1000000L * cluster_id)

m <- arg_int(args, "m", 20L)
p_resp <- as.numeric(arg_value(args, "p_resp", 0.5))
p_cens <- as.numeric(arg_value(args, "p_cens", 0.3))
if (!is.finite(p_resp) || abs(p_resp - 0.5) > 1e-12) {
  stop(
    "The preserved Simulation 2 DGP requires --p-resp 0.5.",
    call. = FALSE
  )
}
if (!is.finite(p_cens) || abs(p_cens - 0.3) > 1e-12) {
  stop(
    "The preserved Simulation 2 DGP requires --p-cens 0.3.",
    call. = FALSE
  )
}
max_runs <- arg_int(
  args,
  "max_runs",
  sim2_profile_default(profile_defaults, "max_runs", 1L)
)
mc_reps <- arg_int(
  args,
  "mc_reps",
  sim2_profile_default(profile_defaults, "mc_reps", 1L)
)
n_trial <- arg_int(args, "n_trial", NA_integer_)
n_aux <- arg_int(args, "n_aux", NA_integer_)

total_cores <- arg_int(args, "total_cores", NA_integer_)
reserve_cores <- arg_int(
  args,
    "reserve_cores",
    sim2_profile_default(profile_defaults, "reserve_cores", 1L)
)
fit_workers <- arg_int(
  args,
  "fit_workers",
  sim2_profile_default(profile_defaults, "fit_workers", NA_integer_)
)
cv_workers_per_fit <- arg_int(
  args,
  "cv_workers_per_fit",
  sim2_profile_default(
    profile_defaults,
    "cv_workers_per_fit",
    NA_integer_
  )
)
max_concurrent_hal_fits <- arg_int(
  args,
  "max_concurrent_hal_fits",
  NA_integer_
)
max_active_mc <- arg_int(args, "max_active_mc", NA_integer_)
active_mc_headroom <- as.numeric(arg_value(
  args,
  "active_mc_headroom",
  1.2
))
cv_folds <- arg_int(
  args,
  "cv_folds",
  sim2_profile_default(profile_defaults, "cv_folds", SIM2_DEFAULT_CV_FOLDS)
)
cv_nlambda <- arg_int(
  args,
  "cv_nlambda",
  sim2_profile_default(profile_defaults, "cv_nlambda", 5L)
)

out_dir <- resolve_cli_path(arg_value(
  args,
  "out_dir",
  file.path(project_dir, "simulation", "sim_data", "sim2")
))
cache_dir <- resolve_cli_path(arg_value(
  args,
  "cache_dir",
  file.path(project_dir, "cache")
))
backend <- arg_value(
  args,
  "backend",
  sim2_profile_default(
    profile_defaults,
    "backend",
    Sys.getenv("HAL_FIT_BACKEND", "auto")
  )
)
use_cache <- sim2_profile_default(profile_defaults, "use_cache", FALSE)
if ("use_cache" %in% names(args)) {
  use_cache <- arg_bool(args, "use_cache", FALSE)
}
if (arg_bool(args, "no_cache", FALSE)) {
  use_cache <- FALSE
}
resume <- !arg_bool(args, "no_resume", FALSE)
checkpoint_dir <- resolve_cli_path(arg_value(args, "checkpoint_dir", NULL))
dry_run <- arg_bool(args, "dry_run", FALSE)
run_analysis <- !arg_bool(args, "no_analysis", FALSE)
figure_dir <- resolve_cli_path(arg_value(
  args,
  "figure_dir",
  file.path(project_dir, "simulation", "sim_figures", "sim2")
))
make_plots <- !arg_bool(args, "no_plots", FALSE)

if (!scale %in% c("demo", "small", "large")) {
  stop("--scale must be one of: demo, small, large.", call. = FALSE)
}
if (!backend %in% c("auto", "fork", "psock")) {
  stop("--backend must be one of: auto, fork, psock.", call. = FALSE)
}
if (!grepl("^[A-Za-z0-9_.-]+$", run_id)) {
  stop(
    "--run-id may contain only letters, digits, underscore, dot, and hyphen.",
    call. = FALSE
  )
}
if (m > 20L || m %% 2L != 0L) {
  stop("--m must be an even integer between 2 and 20.", call. = FALSE)
}

if (is.na(n_trial)) n_trial <- NULL
if (is.na(n_aux)) n_aux <- NULL
if (is.na(max_runs)) max_runs <- NULL
if (is.na(total_cores)) total_cores <- NULL
if (is.na(reserve_cores)) reserve_cores <- NULL
if (is.na(fit_workers)) fit_workers <- NULL
if (is.na(cv_workers_per_fit)) cv_workers_per_fit <- NULL
if (is.na(max_concurrent_hal_fits)) max_concurrent_hal_fits <- NULL
if (is.na(max_active_mc)) max_active_mc <- NULL

if (!is.null(profile_name)) {
  print_sim2_profile_resolution(
    profile_name,
    list(
      scale = scale,
      max_runs = max_runs,
      mc_reps = mc_reps,
      run_id = run_id,
      base_seed = base_seed,
      n_trial = n_trial,
      n_aux = n_aux,
      total_cores = total_cores,
      reserve_cores = reserve_cores,
      fit_workers = fit_workers,
      cv_workers_per_fit = cv_workers_per_fit,
      max_concurrent_hal_fits = max_concurrent_hal_fits,
      max_active_mc = max_active_mc,
      active_mc_headroom = active_mc_headroom,
      cv_folds = cv_folds,
      cv_nlambda = cv_nlambda,
      backend = backend,
      use_cache = use_cache,
      resume = resume
    ),
    explicit_names = names(args)
  )
}


# run simulations ---------------------------------------------------------

run_output <- run_sim2_demo(
  project_dir = project_dir,
  scale = scale,
  run_id = run_id,
  base_seed = base_seed,
  out_dir = out_dir,
  cache_dir = cache_dir,
  use_cache = use_cache,
  resume = resume,
  checkpoint_dir = checkpoint_dir,
  backend = backend,
  m = m,
  p_resp = p_resp,
  p_cens = p_cens,
  n_trial = n_trial,
  n_aux = n_aux,
  max_runs = max_runs,
  mc_reps = mc_reps,
  total_cores = total_cores,
  reserve_cores = reserve_cores,
  fit_workers = fit_workers,
  cv_workers_per_fit = cv_workers_per_fit,
  max_concurrent_hal_fits = max_concurrent_hal_fits,
  max_active_mc = max_active_mc,
  active_mc_headroom = active_mc_headroom,
  cv_folds = cv_folds,
  cv_nlambda = cv_nlambda,
  dry_run = dry_run,
  write_sd_file = TRUE,
  cluster_id = cluster_id
)

if (!isTRUE(dry_run) &&
    (!identical(as.character(run_output$run_summary$status[[1L]]), "ok") ||
     nrow(run_output$timings) != nrow(run_output$sample_grid) ||
     any(as.character(run_output$timings$status) != "ok"))) {
  failed <- run_output$timings$run_index[
    as.character(run_output$timings$status) != "ok"
  ]
  missing <- setdiff(
    seq_len(nrow(run_output$sample_grid)),
    as.integer(run_output$timings$run_index)
  )
  failed <- sort(unique(c(as.integer(failed), missing)))
  stop(
    "One or more simulation runs failed. Failed run_index values: ",
    paste(failed, collapse = ", "),
    ". Timing/checkpoint diagnostics were saved.",
    call. = FALSE
  )
}


# analyze simulations -----------------------------------------------------

if (!isTRUE(dry_run) && isTRUE(run_analysis)) {
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
}
