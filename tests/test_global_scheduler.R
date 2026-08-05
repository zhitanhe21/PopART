###############################################################################
# Fast state-machine and callr integration checks for the global HAL queue
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


# Automatic defaults reproduce the adopted 1.2-headroom rule.
plan4 <- choose_fitcv_plan(
  total_cores = 4L,
  reserve_cores = 0L,
  n_fit_jobs = 4L,
  cv_folds = 3L,
  n_runs = 80L,
  active_mc_headroom = 1.2
)
plan3 <- choose_fitcv_plan(
  total_cores = 4L,
  reserve_cores = 1L,
  n_fit_jobs = 4L,
  cv_folds = 3L,
  n_runs = 80L,
  active_mc_headroom = 1.2
)
stopifnot(
  identical(plan4$fit_workers_per_mc, 2L),
  identical(plan4$global_hal_slots, 4L),
  identical(plan4$max_active_mc, 3L),
  identical(plan3$global_hal_slots, 3L),
  identical(plan3$max_active_mc, 2L)
)

plan2_reserved <- choose_fitcv_plan(
  total_cores = 2L,
  reserve_cores = 1L,
  n_fit_jobs = 4L,
  cv_folds = 3L,
  n_runs = 4L
)
stopifnot(
  identical(plan2_reserved$reserve_cores, 1L),
  identical(plan2_reserved$usable_cores, 1L),
  identical(plan2_reserved$global_hal_slots, 1L)
)


fake_task <- function(spec, use_cache) {
  Sys.sleep(spec$delay)
  if (isTRUE(spec$fail)) stop(spec$error_message)
  list(
    run_id = spec$run_id,
    fit_label = spec$fit_label,
    pid = Sys.getpid()
  )
}

make_fake_prepare <- function(delays, failures = list()) {
  force(delays)
  force(failures)
  function(run_id) {
    labels <- paste0("J", seq_len(4L))
    jobs <- lapply(labels, function(label) {
      failure_key <- paste(run_id, label, sep = ":")
      list(
        run_id = run_id,
        fit_label = label,
        delay = delays[[failure_key]] %||% delays[[as.character(run_id)]] %||% 0.02,
        fail = failure_key %in% names(failures),
        error_message = if (failure_key %in% names(failures)) {
          failures[[failure_key]]
        } else {
          ""
        }
      )
    })
    names(jobs) <- labels
    list(jobs = jobs, payload = run_id)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

fake_finalize <- function(run_id, context, fit_outputs, error, metrics) {
  list(
    run_id = run_id,
    status = if (is.null(error)) "ok" else "error",
    outputs = fit_outputs,
    error = error,
    metrics = metrics
  )
}


# PDF-style laptop case: four global slots, two fits per run, three active runs.
basic <- run_global_hal_scheduler(
  run_ids = seq_len(4L),
  prepare_run = make_fake_prepare(list()),
  finalize_run = fake_finalize,
  global_workers = 4L,
  per_run_workers = 2L,
  max_active_runs = 3L,
  task_fun = fake_task,
  initialize_hal = FALSE,
  verbose = FALSE
)

basic_dispatch <- basic$trace[basic$trace$event == "dispatch", , drop = FALSE]
stopifnot(
  identical(names(basic$finalized), as.character(seq_len(4L))),
  all(vapply(basic$finalized, function(x) identical(x$status, "ok"), logical(1))),
  all(vapply(basic$finalized, function(x) length(x$outputs) == 4L, logical(1))),
  nrow(basic_dispatch) == 16L,
  all(table(basic_dispatch$run_id) == 4L),
  !anyDuplicated(paste(basic_dispatch$run_id, basic_dispatch$fit_label)),
  basic$max_observed_global_fits <= 4L,
  basic$max_observed_active_runs <= 3L,
  basic$max_observed_per_run_fits <= 2L,
  max(basic$trace$running_fits) <= 4L,
  max(basic$trace$active_runs) <= 3L,
  max(basic$trace$max_running_this_run) <= 2L
)

served <- do.call(rbind, lapply(basic$finalized, function(x) {
  do.call(rbind, lapply(x$outputs, function(y) {
    data.frame(run_id = y$run_id, pid = y$pid)
  }))
}))
stopifnot(
  length(unique(served$pid)) <= 4L,
  any(vapply(split(served$run_id, served$pid), function(x) {
    length(unique(x)) > 1L
  }, logical(1)))
)


# Straggler check: run 3 must be admitted while run 1/J1 is still running.
straggler_delays <- list(
  "1:J1" = 0.45,
  "1" = 0.02,
  "2" = 0.02,
  "3" = 0.02
)
straggler <- run_global_hal_scheduler(
  run_ids = 1:3,
  prepare_run = make_fake_prepare(straggler_delays),
  finalize_run = fake_finalize,
  global_workers = 2L,
  per_run_workers = 1L,
  max_active_runs = 2L,
  task_fun = fake_task,
  initialize_hal = FALSE,
  verbose = FALSE
)
trace <- straggler$trace
admit3 <- trace$event_index[trace$event == "admit" & trace$run_id == "3"]
complete_slow <- trace$event_index[
  trace$event == "complete" &
    trace$run_id == "1" &
    trace$fit_label == "J1"
]
stopifnot(
  length(admit3) == 1L,
  length(complete_slow) == 1L,
  admit3 < complete_slow
)


# One ordinary task error fails only its run; other runs continue and finish.
failure <- run_global_hal_scheduler(
  run_ids = 1:3,
  prepare_run = make_fake_prepare(
    list(),
    failures = list("2:J2" = "injected fit failure")
  ),
  finalize_run = fake_finalize,
  global_workers = 3L,
  per_run_workers = 2L,
  max_active_runs = 2L,
  task_fun = fake_task,
  initialize_hal = FALSE,
  verbose = FALSE
)
stopifnot(
  identical(failure$finalized[["1"]]$status, "ok"),
  identical(failure$finalized[["2"]]$status, "error"),
  identical(failure$finalized[["3"]]$status, "ok"),
  grepl("injected fit failure", failure$finalized[["2"]]$error),
  sum(failure$trace$event == "fit_error") == 1L,
  failure$max_observed_global_fits <= 3L,
  failure$max_observed_active_runs <= 2L,
  failure$max_observed_per_run_fits <= 2L
)


# A finalizer exception must be surfaced in the scheduler result so the
# invocation layer cannot silently report a partial aggregate as successful.
finalizer_failure <- run_global_hal_scheduler(
  run_ids = 1L,
  prepare_run = make_fake_prepare(list()),
  finalize_run = function(...) stop("injected finalizer failure"),
  global_workers = 1L,
  per_run_workers = 1L,
  max_active_runs = 1L,
  task_fun = fake_task,
  initialize_hal = FALSE,
  verbose = FALSE
)
stopifnot(
  identical(finalizer_failure$finalized[["1"]]$status, "error"),
  isTRUE(finalizer_failure$finalized[["1"]]$finalizer_failed),
  grepl(
    "injected finalizer failure",
    finalizer_failure$finalized[["1"]]$error
  )
)


# A dead idle slot is retired. Its unsubmitted fit stays pending and the
# surviving persistent slot completes every run instead of poisoning later
# runs that happen to be assigned to the dead node.
dead_slot_pool <- make_global_hal_callr_slots(
  2L,
  initialize_hal = FALSE
)
on.exit(close_global_hal_callr_slots(dead_slot_pool), add = TRUE)
dead_slot_pool[[1L]]$close()
dead_slot_recovery <- run_global_hal_scheduler(
  run_ids = 1:3,
  prepare_run = make_fake_prepare(list()),
  finalize_run = fake_finalize,
  global_workers = 2L,
  per_run_workers = 1L,
  max_active_runs = 2L,
  task_fun = fake_task,
  slots = dead_slot_pool,
  initialize_hal = FALSE,
  verbose = FALSE
)
stopifnot(
  all(vapply(
    dead_slot_recovery$finalized,
    function(x) identical(x$status, "ok"),
    logical(1)
  )),
  sum(dead_slot_recovery$trace$event == "submit_error") == 1L,
  nrow(dead_slot_recovery$trace[
    dead_slot_recovery$trace$event == "dispatch",
    ,
    drop = FALSE
  ]) == 12L
)
close_global_hal_callr_slots(dead_slot_pool)

cat(paste(
  "Global scheduler planning, refill, slot reuse, finalizer,",
  "and failure tests passed.\n"
))
