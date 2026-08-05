###############################################################################
# Installed-package smoke tests (base R only)
###############################################################################

library(aipwLocalDemo)


assert_true <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || !isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}


assert_has_names <- function(object, required, label) {
  missing <- setdiff(required, names(object))
  if (length(missing) > 0L) {
    stop(
      label,
      " is missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# Public API: planning, grids, profiles, data generation, and the dry-run API.
required_exports <- c(
  "choose_fitcv_plan",
  "fitcv_sample_grid",
  "sim2_profile_defaults",
  "simulate_sim2_data",
  "run_sim2_demo"
)
assert_true(
  all(required_exports %in% getNamespaceExports("aipwLocalDemo")),
  "The installed namespace is missing one or more required public functions."
)

plan <- choose_fitcv_plan(
  total_cores = 4L,
  reserve_cores = 1L,
  n_fit_jobs = 4L,
  cv_folds = 3L,
  n_runs = 8L,
  active_mc_headroom = 1.2
)
assert_has_names(
  plan,
  c(
    "total_cores", "reserve_cores", "usable_cores",
    "fit_workers_per_mc", "cv_workers_per_fit",
    "global_hal_slots", "max_active_mc"
  ),
  "CPU plan"
)
stopifnot(
  identical(plan$usable_cores, 3L),
  identical(plan$fit_workers_per_mc, 2L),
  identical(plan$cv_workers_per_fit, 1L),
  identical(plan$global_hal_slots, 3L),
  identical(plan$max_active_mc, 2L)
)

demo_grid <- fitcv_sample_grid("demo")
assert_true(
  identical(demo_grid, fitcv_sample_grid("demo")),
  "fitcv_sample_grid() is not reproducible."
)
stopifnot(
  is.data.frame(demo_grid),
  identical(names(demo_grid), c("n_trial", "n_aux")),
  nrow(demo_grid) == 4L,
  identical(sort(unique(demo_grid$n_trial)), c(120L, 300L)),
  identical(sort(unique(demo_grid$n_aux)), c(120L, 300L))
)

smoke_profile <- sim2_profile_defaults("smoke")
assert_has_names(
  smoke_profile,
  c(
    "scale", "max_runs", "mc_reps", "fit_workers",
    "cv_workers_per_fit", "cv_folds", "cv_nlambda", "use_cache"
  ),
  "Smoke profile"
)
stopifnot(
  identical(smoke_profile$scale, "demo"),
  identical(smoke_profile$max_runs, 1L),
  identical(smoke_profile$mc_reps, 1L),
  identical(smoke_profile$cv_folds, 3L)
)

data_args <- list(
  m = 20L,
  n_trial = 40L,
  n_aux = 40L,
  p_resp = 0.5,
  p_cens = 0.3,
  seed = 20260805L
)
sim_a <- do.call(simulate_sim2_data, data_args)
sim_b <- do.call(simulate_sim2_data, data_args)
assert_true(
  identical(sim_a, sim_b),
  "simulate_sim2_data() did not reproduce an identical seeded data set."
)
assert_has_names(
  sim_a,
  c("dat", "eta0", "eta1", "n_trial", "n_aux"),
  "Simulation output"
)
assert_has_names(
  sim_a$dat,
  c("id", "clust", "X1", "A", "W1", "W2", "W3", "R", "C", "Y", "wt", "S"),
  "Simulated data"
)
stopifnot(
  nrow(sim_a$dat) == 80L,
  identical(as.integer(table(sim_a$dat$S)), c(40L, 40L)),
  is.finite(sim_a$eta0),
  is.finite(sim_a$eta1)
)


# Internal global scheduler: use tiny serialized fake tasks, never initialize
# HAL, and verify that the two-slot ceiling is respected.
run_scheduler <- getFromNamespace(
  "run_global_hal_scheduler",
  "aipwLocalDemo"
)

fake_prepare <- function(run_id) {
  labels <- c("Q0", "Q1")
  jobs <- lapply(labels, function(label) {
    list(run_id = run_id, fit_label = label, delay = 0.01)
  })
  names(jobs) <- labels
  list(jobs = jobs, marker = paste0("run-", run_id))
}

fake_task <- function(spec, use_cache) {
  Sys.sleep(spec$delay)
  list(
    run_id = spec$run_id,
    fit_label = spec$fit_label,
    use_cache = use_cache,
    pid = Sys.getpid()
  )
}

fake_finalize <- function(run_id, context, fit_outputs, error, metrics) {
  list(
    run_id = run_id,
    marker = context$marker,
    status = if (is.null(error)) "ok" else "error",
    outputs = fit_outputs,
    error = error,
    metrics = metrics
  )
}

scheduler_result <- run_scheduler(
  run_ids = 1:3,
  prepare_run = fake_prepare,
  finalize_run = fake_finalize,
  global_workers = 2L,
  per_run_workers = 2L,
  max_active_runs = 2L,
  use_cache = FALSE,
  task_fun = fake_task,
  initialize_hal = FALSE,
  poll_interval_ms = 10L,
  verbose = FALSE
)

dispatch <- scheduler_result$trace[
  scheduler_result$trace$event == "dispatch",
  ,
  drop = FALSE
]
stopifnot(
  length(scheduler_result$finalized) == 3L,
  all(vapply(
    scheduler_result$finalized,
    function(x) identical(x$status, "ok"),
    logical(1)
  )),
  nrow(dispatch) == 6L,
  scheduler_result$max_observed_global_fits <= 2L,
  scheduler_result$max_observed_active_runs <= 2L,
  scheduler_result$max_observed_per_run_fits <= 2L,
  max(scheduler_result$trace$running_fits) <= 2L,
  !anyDuplicated(paste(dispatch$run_id, dispatch$fit_label))
)


# Installed-mode dry run: project_dir=NULL must resolve package resources and
# create only planning output under a disposable directory.
local({
scratch <- tempfile("aipwLocalDemo-package-smoke-")
dir.create(scratch, recursive = TRUE)
on.exit(unlink(scratch, recursive = TRUE, force = TRUE), add = TRUE)

dry <- suppressWarnings(run_sim2_demo(
  project_dir = NULL,
  scale = "demo",
  run_id = "package_smoke",
  out_dir = file.path(scratch, "output"),
  cache_dir = file.path(scratch, "cache"),
  use_cache = FALSE,
  resume = FALSE,
  n_trial = 40L,
  n_aux = 40L,
  max_runs = 1L,
  mc_reps = 2L,
  total_cores = 2L,
  reserve_cores = 0L,
  fit_workers = 1L,
  cv_workers_per_fit = 1L,
  max_concurrent_hal_fits = 1L,
  max_active_mc = 1L,
  cv_folds = 3L,
  cv_nlambda = 5L,
  dry_run = TRUE
))

assert_has_names(dry, c("config", "sample_grid", "files"), "Dry-run result")
assert_has_names(
  dry$config,
  c(
    "usable_cores", "fit_workers_per_mc", "cv_workers_per_fit",
    "global_hal_slots", "max_active_mc", "cv_folds", "cv_nlambda",
    "scheduler", "scheduler_backend"
  ),
  "Dry-run configuration"
)
assert_has_names(
  dry$sample_grid,
  c("n_trial", "n_aux", "mc_rep", "mc_reps", "seed"),
  "Dry-run sample grid"
)
assert_has_names(
  dry$files,
  c(
    "plan", "timing", "fit_timing", "results", "run_summary",
    "scheduler_trace", "checkpoint_dir", "sd"
  ),
  "Dry-run file manifest"
)
stopifnot(
  nrow(dry$sample_grid) == 2L,
  identical(dry$sample_grid$n_trial, rep(40L, 2L)),
  identical(dry$sample_grid$n_aux, rep(40L, 2L)),
  identical(dry$sample_grid$mc_rep, 1:2),
  identical(dry$sample_grid$seed, 91000000L + 1:2),
  identical(dry$config$global_hal_slots, 1L),
  identical(dry$config$max_active_mc, 1L),
  file.exists(dry$files$plan),
  !file.exists(dry$files$results),
  !dir.exists(dry$files$checkpoint_dir)
)

plan_csv <- utils::read.csv(
  dry$files$plan,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
assert_has_names(
  plan_csv,
  c(
    "implementation", "resume", "mc_reps", "total_runs", "seed_design",
    "global_hal_slots", "max_active_mc", "cv_folds", "cv_nlambda",
    "scheduler", "scheduler_backend"
  ),
  "Dry-run plan CSV"
)
stopifnot(
  nrow(plan_csv) == 1L,
  identical(plan_csv$total_runs, 2L),
  identical(plan_csv$mc_reps, 2L)
)
})

cat("Installed-package planning, data, scheduler, and dry-run smoke tests passed.\n")
