###############################################################################
# Lightweight validation for v2 persistence and planning helpers
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
load_fitcvreuse_packages()


plan <- choose_fitcv_plan(
  total_cores = 4L,
  reserve_cores = 2L,
  n_fit_jobs = 4L,
  cv_folds = 4L
)
stopifnot(
  identical(plan$fit_workers, 2L),
  identical(plan$fit_workers_per_mc, 2L),
  identical(plan$cv_workers_per_fit, 1L),
  identical(plan$max_concurrent_hal_fits, 2L),
  identical(plan$max_active_mc, 1L),
  identical(plan$active_mc_headroom, 1.2)
)

mc_grid <- expand_sim2_mc_grid(
  fitcv_sample_grid("small"),
  mc_reps = 20L,
  base_seed = 91000000L
)
combination_counts <- table(mc_grid$n_trial, mc_grid$n_aux)
seed_counts <- table(mc_grid$seed)
stopifnot(
  nrow(mc_grid) == 80L,
  all(combination_counts == 20L),
  length(unique(mc_grid$seed)) == 20L,
  all(seed_counts == 4L),
  identical(unique(mc_grid$seed), 91000001:91000020)
)


config <- list(
  cv_folds = 3L,
  cv_nlambda = 5L,
  cv_workers_per_fit = 1L
)

# Scheduler resource choices must not change the statistical checkpoint
# identity. They affect when and where a fit runs, not what is fitted.
signature_config_a <- utils::modifyList(config, list(
  fit_workers_per_mc = 1L,
  max_concurrent_hal_fits = 1L,
  max_active_mc = 1L
))
signature_config_b <- utils::modifyList(config, list(
  cv_workers_per_fit = 3L,
  fit_workers_per_mc = 2L,
  max_concurrent_hal_fits = 4L,
  max_active_mc = 3L
))
signature_args <- list(
  project_dir = project_dir,
  seed = 93000001L,
  m = 20L,
  n_trial = 120L,
  n_aux = 120L,
  p_resp = 0.5,
  p_cens = 0.3,
  code_fingerprint = "fixed-test-fingerprint"
)
signature_a <- do.call(
  make_run_signature,
  c(signature_args, list(config = signature_config_a))
)
signature_b <- do.call(
  make_run_signature,
  c(signature_args, list(config = signature_config_b))
)
stopifnot(identical(signature_a, signature_b))

sim_a <- simulate_sim2_data(20L, 120L, 120L, 0.5, 0.3, 93000001L)
sim_b <- simulate_sim2_data(20L, 120L, 120L, 0.5, 0.3, 93000002L)

temp_root <- tempfile("aipw_local_demo_2_test_")
dir.create(temp_root, recursive = TRUE)
on.exit(unlink(temp_root, recursive = TRUE), add = TRUE)

jobs_a <- make_fitcvreuse_jobs(
  sim_a$dat,
  93000001L,
  sim_a$n_trial,
  sim_a$n_aux,
  config,
  file.path(temp_root, "cache_a")
)
jobs_b <- make_fitcvreuse_jobs(
  sim_b$dat,
  93000002L,
  sim_b$n_trial,
  sim_b$n_aux,
  config,
  file.path(temp_root, "cache_b")
)
config_parallel_cv <- config
config_parallel_cv$cv_workers_per_fit <- 2L
jobs_a_parallel_cv <- make_fitcvreuse_jobs(
  sim_a$dat,
  93000001L,
  sim_a$n_trial,
  sim_a$n_aux,
  config_parallel_cv,
  file.path(temp_root, "cache_a_parallel_cv")
)
keys_a <- vapply(jobs_a$jobs, function(x) x$cache_key, character(1))
keys_b <- vapply(jobs_b$jobs, function(x) x$cache_key, character(1))
keys_a_parallel_cv <- vapply(
  jobs_a_parallel_cv$jobs,
  function(x) x$cache_key,
  character(1)
)
stopifnot(
  !identical(keys_a, keys_b),
  identical(keys_a, keys_a_parallel_cv),
  all(vapply(
    jobs_a$jobs,
    function(x) identical(x$num_knots, c(5, 3)),
    logical(1)
  )),
  all(vapply(
    jobs_a$jobs,
    function(x) identical(x$cv_folds, 3L),
    logical(1)
  ))
)


rds_path <- file.path(temp_root, "atomic", "object.rds")
csv_path <- file.path(temp_root, "atomic", "table.csv")
atomic_save_rds(list(value = 42L), rds_path)
atomic_write_csv(data.frame(value = 42L), csv_path)
stopifnot(
  identical(readRDS(rds_path)$value, 42L),
  identical(utils::read.csv(csv_path)$value, 42L)
)

invalid_cache_path <- file.path(temp_root, "atomic", "invalid_cache.rds")
saveRDS(list(not_a_fit = TRUE), invalid_cache_path)
suppressWarnings(
  invalid_cache <- read_fit_cache(invalid_cache_path, "expected-key")
)
stopifnot(
  is.null(invalid_cache),
  !file.exists(invalid_cache_path),
  length(Sys.glob(paste0(invalid_cache_path, ".invalid-fit-cache-*"))) == 1L
)


signature <- hash_r_object(list(test = "checkpoint"))
checkpoint_dir <- file.path(temp_root, "checkpoints")
attempt_id <- "test-attempt"
checkpoint_path <- make_run_checkpoint_file(
  checkpoint_dir,
  1L,
  signature,
  attempt_id
)
checkpoint <- list(
  schema = RUN_CHECKPOINT_SCHEMA,
  run_signature = signature,
  attempt_id = attempt_id,
  status = "ok",
  timing = data.frame(status = "ok"),
  fit_timings = data.frame(fit_label = EXPECTED_HAL_FIT_LABELS),
  results = data.frame(
    Estimator = rep("AIPW", 2L),
    Version = c("Naive", "Proposed")
  )
)
atomic_save_rds(checkpoint, checkpoint_path)
restored <- read_successful_run_checkpoint(checkpoint_path, signature)
found <- find_successful_run_checkpoint(checkpoint_dir, 1L, signature)
stopifnot(
  !is.null(restored),
  identical(found$path, checkpoint_path),
  nrow(restored$fit_timings) == 4L,
  nrow(restored$results) == 2L
)

invisible(write_run_output_shards(
  restored,
  checkpoint_dir,
  1L,
  signature,
  attempt_id
))
shards <- list.files(
  file.path(checkpoint_dir, "csv"),
  pattern = "[.]csv$",
  full.names = TRUE
)
stopifnot(length(shards) == 3L)


recovery_path <- file.path(temp_root, "atomic", "recover.rds")
recovery_backup <- paste0(recovery_path, ".previous-test")
saveRDS(list(value = "recoverable"), recovery_backup)
suppressWarnings(recovered <- recover_atomic_target(recovery_path))
stopifnot(
  isTRUE(recovered),
  identical(readRDS(recovery_path)$value, "recoverable")
)


fingerprint_project <- file.path(temp_root, "fingerprint_project")
dir.create(fingerprint_project, recursive = TRUE)
invisible(file.copy(
  file.path(project_dir, "R"),
  fingerprint_project,
  recursive = TRUE
))
dir.create(file.path(fingerprint_project, "tests"))
fingerprint_before <- make_code_fingerprint(fingerprint_project)
irrelevant_test_file <- file.path(
  fingerprint_project,
  "tests",
  "fingerprint_probe.R"
)
writeLines("# presentation-only fingerprint probe", irrelevant_test_file)
fingerprint_after <- make_code_fingerprint(fingerprint_project)
stopifnot(identical(fingerprint_before, fingerprint_after))


cat("All persistence and planning helper checks passed.\n")
