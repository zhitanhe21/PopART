###############################################################################
# Lightweight validation for Simulation 2 convenience profiles
###############################################################################


args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) > 0L) {
  test_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  test_dir <- normalizePath("tests")
}
project_dir <- normalizePath(file.path(test_dir, ".."))

source(file.path(project_dir, "R", "source_all.R"))
source_aipw_local_demo(project_dir)


assert_profile <- function(name, expected) {
  actual <- sim2_profile_defaults(name)
  stopifnot(setequal(names(actual), names(expected)))
  for (field in names(expected)) {
    stopifnot(identical(actual[[field]], expected[[field]]))
  }
}


assert_profile(
  "smoke",
  list(
    scale = "demo",
    max_runs = 1L,
    mc_reps = 1L,
    reserve_cores = 0L,
    fit_workers = 1L,
    cv_workers_per_fit = 1L,
    cv_folds = 3L,
    cv_nlambda = 50L,
    backend = "auto",
    use_cache = FALSE
  )
)

common_small <- list(
  scale = "small",
  max_runs = 4L,
  reserve_cores = 0L,
  fit_workers = 2L,
  cv_workers_per_fit = 1L,
  cv_folds = 3L,
  cv_nlambda = 50L,
  backend = "psock",
  use_cache = FALSE
)
assert_profile(
  "quick",
  utils::modifyList(common_small, list(mc_reps = 1L), keep.null = TRUE)
)
assert_profile(
  "hour",
  utils::modifyList(common_small, list(mc_reps = 2L), keep.null = TRUE)
)
assert_profile(
  "formal",
  utils::modifyList(common_small, list(mc_reps = 20L), keep.null = TRUE)
)

unknown_profile <- try(
  sim2_profile_defaults("overnight"),
  silent = TRUE
)
stopifnot(inherits(unknown_profile, "try-error"))


resolved <- resolve_sim2_profile_args(
  "formal",
  overrides = list(
    max_runs = 2L,
    mc_reps = 3L,
    fit_workers = 2L,
    cv_nlambda = 7L,
    use_cache = TRUE
  ),
  legacy_defaults = list(
    scale = "demo",
    max_runs = 1L,
    mc_reps = 1L,
    reserve_cores = 2L,
    cv_nlambda = 5L,
    use_cache = FALSE
  )
)
stopifnot(
  identical(resolved$scale, "small"),
  identical(resolved$max_runs, 2L),
  identical(resolved$mc_reps, 3L),
  identical(resolved$reserve_cores, 0L),
  identical(resolved$fit_workers, 2L),
  identical(resolved$cv_nlambda, 7L),
  isTRUE(resolved$use_cache)
)

explicit_null <- resolve_sim2_profile_args(
  "smoke",
  overrides = list(max_runs = NULL)
)
stopifnot(
  "max_runs" %in% names(explicit_null),
  is.null(explicit_null$max_runs)
)


temp_root <- tempfile("aipw_local_demo_2_profiles_")
dir.create(temp_root, recursive = TRUE)
on.exit(unlink(temp_root, recursive = TRUE), add = TRUE)

run_profile_dry <- function(profile) {
  output <- NULL
  invisible(capture.output(
    output <- suppressWarnings(run_sim2_demo_profile(
      project_dir = project_dir,
      profile = profile,
      run_id = paste0("test_", profile),
      out_dir = file.path(temp_root, profile),
      dry_run = TRUE
    ))
  ))
  output
}

profile_outputs <- lapply(
  SIM2_PROFILE_NAMES,
  run_profile_dry
)
names(profile_outputs) <- SIM2_PROFILE_NAMES

expected_rows <- c(
  smoke = 1L,
  quick = 4L,
  hour = 8L,
  formal = 80L
)
for (profile in SIM2_PROFILE_NAMES) {
  output <- profile_outputs[[profile]]
  defaults <- sim2_profile_defaults(profile)
  stopifnot(
    identical(nrow(output$sample_grid), expected_rows[[profile]]),
    identical(output$profile, profile),
    file.exists(output$files$plan),
    !dir.exists(output$files$checkpoint_dir)
  )
  for (field in names(defaults)) {
    stopifnot(
      identical(output$profile_arguments[[field]], defaults[[field]])
    )
  }
  expected_slots <- min(
    floor(output$config$usable_cores / output$config$cv_workers_per_fit),
    nrow(output$sample_grid) * output$config$fit_workers_per_mc
  )
  expected_active <- min(
    nrow(output$sample_grid),
    ceiling(
      1.2 * expected_slots / output$config$fit_workers_per_mc
    )
  )
  stopifnot(
    identical(output$config$scheduler, "global_hal_queue"),
    identical(output$config$scheduler_backend, "callr"),
    identical(output$config$max_concurrent_hal_fits, as.integer(expected_slots)),
    identical(output$config$max_active_mc, as.integer(expected_active))
  )
}

formal_grid <- profile_outputs$formal$sample_grid
formal_combo_counts <- table(formal_grid$n_trial, formal_grid$n_aux)
formal_seed_counts <- table(formal_grid$seed)
stopifnot(
  all(formal_combo_counts == 20L),
  length(unique(formal_grid$seed)) == 20L,
  all(formal_seed_counts == 4L),
  identical(unique(formal_grid$seed), 91000001:91000020)
)


override_output <- NULL
invisible(capture.output(
  override_output <- suppressWarnings(run_sim2_demo_profile(
    project_dir = project_dir,
    profile = "formal",
    max_runs = 2L,
    mc_reps = 3L,
    fit_workers = 2L,
    cv_nlambda = 7L,
    use_cache = TRUE,
    run_id = "test_formal_override",
    out_dir = file.path(temp_root, "formal_override"),
    cache_dir = file.path(temp_root, "formal_override_cache"),
    dry_run = TRUE
  ))
))
override_seed_counts <- table(override_output$sample_grid$seed)
stopifnot(
  nrow(override_output$sample_grid) == 6L,
  all(override_seed_counts == 2L),
  override_output$config$fit_workers <= 2L,
  identical(override_output$config$cv_nlambda, 7L),
  isTRUE(override_output$profile_arguments$use_cache),
  !dir.exists(override_output$files$checkpoint_dir)
)


quick_explicit <- sim2_profile_defaults("quick")
explicit_output <- NULL
invisible(capture.output(
  explicit_output <- suppressWarnings(run_sim2_demo(
    project_dir = project_dir,
    scale = quick_explicit$scale,
    run_id = "test_quick_explicit",
    out_dir = file.path(temp_root, "quick_explicit"),
    use_cache = quick_explicit$use_cache,
    backend = quick_explicit$backend,
    max_runs = quick_explicit$max_runs,
    mc_reps = quick_explicit$mc_reps,
    reserve_cores = quick_explicit$reserve_cores,
    fit_workers = quick_explicit$fit_workers,
    cv_workers_per_fit = quick_explicit$cv_workers_per_fit,
    cv_folds = quick_explicit$cv_folds,
    cv_nlambda = quick_explicit$cv_nlambda,
    dry_run = TRUE
  ))
))

quick_profile <- profile_outputs$quick
config_fields <- c(
  "total_cores",
  "reserve_cores",
  "usable_cores",
  "fit_workers",
  "fit_workers_per_mc",
  "cv_workers_per_fit",
  "max_concurrent_hal_fits",
  "max_active_mc",
  "active_mc_headroom",
  "cv_folds",
  "cv_nlambda",
  "backend",
  "scheduler",
  "scheduler_backend"
)
stopifnot(
  identical(quick_profile$sample_grid, explicit_output$sample_grid),
  identical(
    quick_profile$config[config_fields],
    explicit_output$config[config_fields]
  )
)

code_fingerprint <- make_code_fingerprint(project_dir)
profile_signature <- make_run_signature(
  project_dir = project_dir,
  seed = quick_profile$sample_grid$seed[[1]],
  m = 20L,
  n_trial = quick_profile$sample_grid$n_trial[[1]],
  n_aux = quick_profile$sample_grid$n_aux[[1]],
  p_resp = 0.5,
  p_cens = 0.3,
  config = quick_profile$config,
  code_fingerprint = code_fingerprint
)
explicit_signature <- make_run_signature(
  project_dir = project_dir,
  seed = explicit_output$sample_grid$seed[[1]],
  m = 20L,
  n_trial = explicit_output$sample_grid$n_trial[[1]],
  n_aux = explicit_output$sample_grid$n_aux[[1]],
  p_resp = 0.5,
  p_cens = 0.3,
  config = explicit_output$config,
  code_fingerprint = code_fingerprint
)
stopifnot(identical(profile_signature, explicit_signature))

cat("Profile helper and dry-run tests passed.\n")
