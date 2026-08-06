###############################################################################
###############################################################################

# PopART manuscript-reference Monte Carlo workflow

###############################################################################
###############################################################################


#' Set an integer environment variable when a value is provided.
#'
#' @param value Value to coerce to integer. `NULL` or `NA` values are ignored.
#' @param env_name Character environment-variable name.
#'
#' @return Invisibly returns `TRUE` if the variable was set and `FALSE`
#'   otherwise.
set_int_env_value <- function(value, env_name) {
  if (is.null(value)) {
    return(invisible(FALSE))
  }

  value <- suppressWarnings(as.integer(value))
  if (is.na(value)) {
    return(invisible(FALSE))
  }

  do.call(Sys.setenv, as.list(stats::setNames(as.character(value), env_name)))
  invisible(TRUE)
}


#' Set a positive numeric environment variable when a value is provided.
#' @noRd
set_numeric_env_value <- function(value, env_name) {
  if (is.null(value)) return(invisible(FALSE))
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    return(invisible(FALSE))
  }
  do.call(Sys.setenv, as.list(stats::setNames(as.character(value), env_name)))
  invisible(TRUE)
}


#' Refresh invocation metadata on rows restored from a checkpoint.
#' @noRd
refresh_resumed_metadata <- function(data, run_id, sample_design, run_index,
                                     total_runs, mc_rep, mc_reps, cpu_plan,
                                     config, use_cache, checkpoint_file) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    return(data)
  }

  source_fields <- c(
    "run_id", "sample_design", "total_runs", "mc_rep", "mc_reps",
    "cpu_plan", "total_cores",
    "reserve_cores", "usable_cores", "fit_workers",
    "fit_workers_per_mc", "cv_workers_per_fit",
    "max_concurrent_hal_fits", "max_active_mc", "active_mc_headroom",
    "scheduler",
    "scheduler_backend", "backend", "legacy_fit_backend", "use_cache",
    "checkpoint_file"
  )
  for (field in source_fields) {
    if (field %in% names(data)) {
      data[[paste0("checkpoint_source_", field)]] <- data[[field]]
    }
  }

  current <- list(
    run_id = run_id,
    sample_design = sample_design,
    implementation = REFERENCE_STUDY_IMPLEMENTATION,
    cpu_plan = cpu_plan,
    run_index = as.integer(run_index),
    total_runs = as.integer(total_runs),
    mc_rep = as.integer(mc_rep),
    mc_reps = as.integer(mc_reps),
    total_cores = config$total_cores,
    reserve_cores = config$reserve_cores,
    usable_cores = config$usable_cores,
    fit_workers = config$fit_workers,
    fit_workers_per_mc = config$fit_workers_per_mc,
    cv_workers_per_fit = config$cv_workers_per_fit,
    max_concurrent_hal_fits = config$max_concurrent_hal_fits,
    max_active_mc = config$max_active_mc,
    active_mc_headroom = config$active_mc_headroom,
    scheduler = config$scheduler,
    scheduler_backend = config$scheduler_backend,
    backend = config$scheduler_backend,
    legacy_fit_backend = config$backend,
    use_cache = isTRUE(use_cache),
    resumed = TRUE,
    checkpoint_file = checkpoint_file
  )
  for (field in names(current)) {
    data[[field]] <- current[[field]]
  }
  data
}


#' Run the PopART manuscript-reference Monte Carlo study.
#'
#' This repository-only entry point uses a bounded cross-replicate HAL queue,
#' worker slots, content-addressed fit caching, successful-run resume
#' checkpoints, and per-run atomic CSV shards.
#'
#' @param sample_sizes Nonempty data frame with explicit `n_trial` and
#'   `n_auxiliary` columns. The study never infers sample sizes from labels such
#'   as small, medium, or large.
#' @param project_dir Optional character path to the repository root.
#' @param run_id Character identifier used in output file names.
#' @param base_seed Integer base seed. Monte Carlo replicate `r` uses
#'   `base_seed + r` across all requested sample-size combinations, matching the
#'   manuscript reference-study grid organization.
#' @param mc_reps Positive integer number of Monte Carlo replicates for every
#'   requested sample-size combination.
#' @param out_dir Character directory for result and timing files.
#' @param cache_dir Character directory for optional HAL fit cache files.
#' @param use_cache Logical, whether HAL fit caching is enabled.
#' @param resume Logical, whether valid successful-run checkpoints may be used.
#' @param checkpoint_dir Optional checkpoint directory. The default is a
#'   run-specific directory below `out_dir/checkpoints`.
#' @param backend Character fit-level backend: `auto`, `fork`, or `psock`.
#' @param m Integer number of clusters.
#' @param p_resp Numeric response probability setting, fixed at `0.5` by the
#'   manuscript reference design.
#' @param p_cens Numeric censoring probability setting, fixed at `0.3` by the
#'   manuscript reference design.
#' @param max_runs Optional integer limiting how many sample-size combinations
#'   are selected before Monte Carlo expansion.
#' @param total_cores,reserve_cores,fit_workers,cv_workers_per_fit,cv_folds,cv_nlambda
#'   Optional local parallel and HAL CV overrides. In the global scheduler,
#'   `fit_workers` is the per-run HAL-fit cap.
#' @param max_concurrent_hal_fits Optional global complete-HAL-fit slot limit.
#' @param max_active_mc Optional bound on prepared unfinished run contexts.
#' @param active_mc_headroom Optional positive multiplier for the automatic
#'   active-run window; the default is `1.2`.
#' @param dry_run Logical, whether to print the plan without fitting HAL models.
#' @return Invisibly returns configuration, paths, timings, and results.
#'
run_monte_carlo_study <- function(
    sample_sizes,
    project_dir = NULL,
    run_id = paste0("study_", format(Sys.time(), "%Y%m%d_%H%M%S")),
    base_seed = 91000000L,
    out_dir = file.path(getwd(), "simulation", "results"),
    cache_dir = file.path(
      tools::R_user_dir("popart", which = "cache"),
      "hal"
    ),
    use_cache = TRUE,
    resume = TRUE,
    checkpoint_dir = NULL,
    backend = Sys.getenv("HAL_FIT_BACKEND", "auto"),
    m = 20L,
    p_resp = 0.5,
    p_cens = 0.3,
    max_runs = NULL,
    mc_reps = 1L,
    total_cores = NULL,
    reserve_cores = NULL,
    fit_workers = NULL,
    cv_workers_per_fit = NULL,
    max_concurrent_hal_fits = NULL,
    max_active_mc = NULL,
    active_mc_headroom = NULL,
    cv_folds = NULL,
    cv_nlambda = NULL,
    dry_run = FALSE) {

  invocation_start_time <- Sys.time()
  if (is.null(project_dir) || !nzchar(project_dir)) {
    project_dir <- getwd()
  }
  project_dir <- normalizePath(project_dir)
  if (length(run_id) != 1L || !nzchar(run_id) ||
      !grepl("^[A-Za-z0-9_.-]+$", run_id)) {
    stop(
      "run_id may contain only letters, digits, underscore, dot, and hyphen.",
      call. = FALSE
    )
  }
  m <- suppressWarnings(as.integer(m))
  if (is.na(m) || m < 2L || m > 20L || m %% 2L != 0L) {
    stop("m must be an even integer between 2 and 20.", call. = FALSE)
  }
  base_seed <- suppressWarnings(as.integer(base_seed))
  if (is.na(base_seed) || base_seed < 0L) {
    stop("base_seed must be a nonnegative integer.", call. = FALSE)
  }
  mc_reps <- suppressWarnings(as.integer(mc_reps))
  if (is.na(mc_reps) || mc_reps < 1L) {
    stop("mc_reps must be a positive integer.", call. = FALSE)
  }
  p_resp <- suppressWarnings(as.numeric(p_resp))
  p_cens <- suppressWarnings(as.numeric(p_cens))
  if (length(p_resp) != 1L || !is.finite(p_resp) ||
      abs(p_resp - 0.5) > 1e-12) {
    stop(
      "The manuscript reference design requires p_resp = 0.5.",
      call. = FALSE
    )
  }
  if (length(p_cens) != 1L || !is.finite(p_cens) ||
      abs(p_cens - 0.3) > 1e-12) {
    stop(
      "The manuscript reference design requires p_cens = 0.3.",
      call. = FALSE
    )
  }
  hal_env_names <- c(
    "HAL_TOTAL_CORES",
    "HAL_RESERVE_CORES",
    "HAL_FIT_WORKERS",
    "HAL_CV_WORKERS_PER_FIT",
    "HAL_MAX_CONCURRENT_HAL_FITS",
    "HAL_MAX_ACTIVE_MC",
    "HAL_ACTIVE_MC_HEADROOM",
    "HAL_CV_FOLDS",
    "HAL_CV_NLAMBDA"
  )
  math_env_names <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  managed_env_names <- c(hal_env_names, math_env_names)
  old_hal_env <- vapply(
    managed_env_names,
    Sys.getenv,
    character(1),
    unset = NA_character_
  )
  on.exit({
    unset_names <- names(old_hal_env)[is.na(old_hal_env)]
    if (length(unset_names) > 0L) {
      Sys.unsetenv(unset_names)
    }
    set_values <- old_hal_env[!is.na(old_hal_env)]
    if (length(set_values) > 0L) {
      do.call(Sys.setenv, as.list(set_values))
    }
  }, add = TRUE)

  manual_worker_plan <- !is.null(fit_workers) ||
    !is.null(cv_workers_per_fit) ||
    !is.null(max_concurrent_hal_fits) ||
    !is.null(max_active_mc) ||
    nzchar(Sys.getenv("HAL_FIT_WORKERS", unset = "")) ||
    nzchar(Sys.getenv("HAL_CV_WORKERS_PER_FIT", unset = "")) ||
    nzchar(Sys.getenv("HAL_MAX_CONCURRENT_HAL_FITS", unset = "")) ||
    nzchar(Sys.getenv("HAL_MAX_ACTIVE_MC", unset = ""))
  cpu_plan_label <- if (isTRUE(manual_worker_plan)) {
    "globalmanual_v2"
  } else {
    "globalauto_v2"
  }

  set_int_env_value(total_cores, "HAL_TOTAL_CORES")
  set_int_env_value(reserve_cores, "HAL_RESERVE_CORES")
  set_int_env_value(fit_workers, "HAL_FIT_WORKERS")
  set_int_env_value(cv_workers_per_fit, "HAL_CV_WORKERS_PER_FIT")
  set_int_env_value(
    max_concurrent_hal_fits,
    "HAL_MAX_CONCURRENT_HAL_FITS"
  )
  set_int_env_value(max_active_mc, "HAL_MAX_ACTIVE_MC")
  set_numeric_env_value(active_mc_headroom, "HAL_ACTIVE_MC_HEADROOM")
  set_int_env_value(cv_folds, "HAL_CV_FOLDS")
  set_int_env_value(cv_nlambda, "HAL_CV_NLAMBDA")

  set_single_threaded_math()
  load_reference_study_packages()

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(use_cache)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    cache_dir <- NULL
  }

  if (!is.data.frame(sample_sizes) || nrow(sample_sizes) < 1L ||
      !all(c("n_trial", "n_auxiliary") %in% names(sample_sizes))) {
    stop(
      "sample_sizes must be a nonempty data frame with n_trial and n_auxiliary columns.",
      call. = FALSE
    )
  }
  sizes <- data.frame(
    n_trial = as.integer(sample_sizes$n_trial),
    n_aux = as.integer(sample_sizes$n_auxiliary)
  )
  sample_design <- "user_supplied"

  if (!is.null(max_runs) && !is.na(as.integer(max_runs))) {
    max_runs <- as.integer(max_runs)
    if (max_runs < 1L) {
      stop("max_runs must be at least 1.", call. = FALSE)
    }
    sizes <- sizes[
      seq_len(min(max_runs, nrow(sizes))),
      ,
      drop = FALSE
    ]
  }
  if (any(is.na(sizes$n_trial)) || any(is.na(sizes$n_aux)) ||
      any(sizes$n_trial < m) || any(sizes$n_aux < m)) {
    stop("Every n_trial and n_aux value must be at least m.", call. = FALSE)
  }
  sizes <- expand_monte_carlo_grid(
    sizes,
    mc_reps = mc_reps,
    base_seed = base_seed
  )

  config <- configure_fitcv_parallel(
    n_fit_jobs = 4L,
    default_cv_folds = REFERENCE_DEFAULT_CV_FOLDS,
    default_cv_nlambda = 50L,
    backend = backend,
    n_runs = nrow(sizes),
    default_active_mc_headroom = 1.2
  )
  config$backend <- resolve_fit_backend(config$backend)
  config$scheduler <- "global_hal_queue"
  config$scheduler_backend <- "callr"

  prefix <- paste("monte_carlo", run_id, sep = "_")
  timing_file <- file.path(out_dir, paste0(prefix, "_timing.csv"))
  fit_timing_file <- file.path(out_dir, paste0(prefix, "_fit_timing.csv"))
  result_file <- file.path(out_dir, paste0(prefix, "_results.csv"))
  plan_file <- file.path(out_dir, paste0(prefix, "_parallel_plan.csv"))
  summary_file <- file.path(out_dir, paste0(prefix, "_run_summary.csv"))
  scheduler_trace_file <- file.path(
    out_dir,
    paste0(prefix, "_scheduler_trace.csv")
  )
  if (is.null(checkpoint_dir) || !nzchar(checkpoint_dir)) {
    checkpoint_dir <- file.path(out_dir, "checkpoints", prefix)
  }

  plan_output <- as.data.frame(config, stringsAsFactors = FALSE)
  plan_output$implementation <- REFERENCE_STUDY_IMPLEMENTATION
  plan_output$resume <- isTRUE(resume)
  plan_output$mc_reps <- mc_reps
  plan_output$sample_size_combinations <- length(unique(paste(
    sizes$n_trial,
    sizes$n_aux,
    sep = "/"
  )))
  plan_output$total_runs <- nrow(sizes)
  plan_output$seed_design <- "shared_by_mc_rep_across_sample_sizes"
  atomic_write_csv(plan_output, plan_file, overwrite = TRUE)

  cat("Running the PopART manuscript-reference Monte Carlo study\n")
  cat("  run_id = ", run_id, "\n", sep = "")
  cat("  base_seed = ", base_seed, "\n", sep = "")
  cat("  mc_reps = ", mc_reps, "\n", sep = "")
  cat("  total_runs = ", nrow(sizes), "\n", sep = "")
  cat("  m = ", m, "\n", sep = "")
  cat("  p_resp = ", p_resp, "\n", sep = "")
  cat("  p_cens = ", p_cens, "\n", sep = "")
  cat("  output_dir = ", out_dir, "\n", sep = "")
  cat(
    "  cache_dir = ",
    if (is.null(cache_dir)) "disabled" else cache_dir,
    "\n",
    sep = ""
  )
  cat("  resume = ", isTRUE(resume), "\n", sep = "")
  cat("  checkpoint_dir = ", checkpoint_dir, "\n", sep = "")
  cat("  timing_file = ", timing_file, "\n", sep = "")
  cat("  fit_timing_file = ", fit_timing_file, "\n", sep = "")
  cat("  result_file = ", result_file, "\n", sep = "")
  cat("  scheduler_trace_file = ", scheduler_trace_file, "\n", sep = "")
  print_parallel_plan(config)
  flush.console()

  output_files <- list(
    plan = plan_file,
    timing = timing_file,
    fit_timing = fit_timing_file,
    results = result_file,
    run_summary = summary_file,
    scheduler_trace = scheduler_trace_file,
    checkpoint_dir = checkpoint_dir
  )

  if (isTRUE(dry_run)) {
    cat("Planned sample grid\n")
    print(sizes)
    cat("Dry run requested; no worker pool or HAL fit was started.\n")
    return(invisible(list(
      config = config,
      sample_grid = sizes,
      files = output_files
    )))
  }

  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  code_fingerprint <- make_code_fingerprint(project_dir)

  timings <- vector("list", nrow(sizes))
  fit_timings_all <- vector("list", nrow(sizes))
  results <- vector("list", nrow(sizes))

  scheduler_trace <- data.frame()
  scheduler_stats <- list(
    max_observed_global_fits = 0L,
    max_observed_active_runs = 0L,
    max_observed_per_run_fits = 0L
  )

  if (identical(config$scheduler, "global_hal_queue")) {
    run_signatures <- character(nrow(sizes))
    pending_indices <- integer()
    resume_trace_rows <- list()

    for (ii in seq_len(nrow(sizes))) {
      seed <- sizes$seed[[ii]]
      mc_rep_ii <- sizes$mc_rep[[ii]]
      run_signature <- make_run_signature(
        project_dir = project_dir,
        seed = seed,
        m = m,
        n_trial = sizes$n_trial[[ii]],
        n_aux = sizes$n_aux[[ii]],
        p_resp = p_resp,
        p_cens = p_cens,
        config = config,
        code_fingerprint = code_fingerprint
      )
      run_signatures[[ii]] <- run_signature
      checkpoint_match <- if (isTRUE(resume)) {
        find_successful_run_checkpoint(
          checkpoint_dir,
          ii,
          run_signature
        )
      } else {
        list(checkpoint = NULL, path = NULL)
      }
      checkpoint <- checkpoint_match$checkpoint

      if (is.null(checkpoint)) {
        pending_indices <- c(pending_indices, ii)
        next
      }

      checkpoint_file <- checkpoint_match$path
      write_run_output_shards(
        checkpoint,
        checkpoint_dir,
        ii,
        run_signature,
        attempt_id = checkpoint$attempt_id
      )
      checkpoint$timing <- refresh_resumed_metadata(
        checkpoint$timing,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )
      checkpoint$fit_timings <- refresh_resumed_metadata(
        checkpoint$fit_timings,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )
      checkpoint$results <- refresh_resumed_metadata(
        checkpoint$results,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )
      timings[[ii]] <- checkpoint$timing
      fit_timings_all[[ii]] <- checkpoint$fit_timings
      results[[ii]] <- checkpoint$results
      resume_trace_rows[[length(resume_trace_rows) + 1L]] <- data.frame(
        event_index = length(resume_trace_rows) + 1L,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"),
        elapsed_sec = as.numeric(difftime(
          Sys.time(), invocation_start_time, units = "secs"
        )),
        event = "resume",
        run_id = as.character(ii),
        fit_label = NA_character_,
        node = NA_integer_,
        task_id = NA_character_,
        active_runs = 0L,
        running_fits = 0L,
        free_slots = NA_integer_,
        max_running_this_run = 0L,
        message = paste("restored from", basename(checkpoint_file)),
        stringsAsFactors = FALSE
      )
      cat(sprintf(
        "[%s] resume ii=%d/%d mc=%d/%d seed=%d from %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        seed,
        basename(checkpoint_file)
      ))
      flush.console()
    }

    prepare_global_run <- function(ii) {
      seed <- sizes$seed[[ii]]
      cat(sprintf(
        "[%s] prepare ii=%d/%d mc=%d/%d seed=%d n_trial=%d n_aux=%d\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        ii,
        nrow(sizes),
        sizes$mc_rep[[ii]],
        mc_reps,
        seed,
        sizes$n_trial[[ii]],
        sizes$n_aux[[ii]]
      ))
      flush.console()
      prepare_reference_replicate(
        m = m,
        n_trial = sizes$n_trial[[ii]],
        n_aux = sizes$n_aux[[ii]],
        p_resp = p_resp,
        p_cens = p_cens,
        seed = seed,
        config = config,
        cache_dir = cache_dir
      )
    }

    finalize_global_run <- function(ii, context, fit_outputs, error, metrics) {
      seed <- sizes$seed[[ii]]
      mc_rep_ii <- sizes$mc_rep[[ii]]
      n_trial_ii <- sizes$n_trial[[ii]]
      n_aux_ii <- sizes$n_aux[[ii]]
      status <- if (is.null(error)) "ok" else "error"
      err_msg <- if (is.null(error)) NA_character_ else as.character(error)
      assembled <- NULL

      if (identical(status, "ok")) {
        assembled <- tryCatch(
          assemble_reference_replicate(context, fit_outputs),
          error = function(e) e
        )
        if (inherits(assembled, "error")) {
          status <- "error"
          err_msg <- paste("AIPW assembly failed:", conditionMessage(assembled))
          assembled <- NULL
        }
      }

      elapsed <- as.numeric(difftime(
        Sys.time(),
        metrics$admitted_at,
        units = "secs"
      ))
      attempt_id <- make_run_attempt_id()
      run_signature <- run_signatures[[ii]]
      checkpoint_file <- make_run_checkpoint_file(
        checkpoint_dir,
        ii,
        run_signature,
        attempt_id = attempt_id
      )

      timing_row <- data.frame(
        run_id = run_id,
        sample_design = sample_design,
        implementation = REFERENCE_STUDY_IMPLEMENTATION,
        cpu_plan = cpu_plan_label,
        scheduler = config$scheduler,
        scheduler_backend = config$scheduler_backend,
        run_index = ii,
        total_runs = nrow(sizes),
        mc_rep = mc_rep_ii,
        mc_reps = mc_reps,
        seed = seed,
        n_trial = n_trial_ii,
        n_aux = n_aux_ii,
        status = status,
        elapsed_sec = elapsed,
        elapsed_min = elapsed / 60,
        prepare_elapsed_sec = metrics$prepare_elapsed_sec,
        queue_wait_sec = metrics$queue_wait_sec,
        fit_span_sec = metrics$fit_span_sec,
        result_rows = if (is.null(assembled)) {
          NA_integer_
        } else {
          nrow(assembled$results)
        },
        error = err_msg,
        total_cores = config$total_cores,
        reserve_cores = config$reserve_cores,
        usable_cores = config$usable_cores,
        fit_workers = config$fit_workers,
        fit_workers_per_mc = config$fit_workers_per_mc,
        max_concurrent_hal_fits = config$max_concurrent_hal_fits,
        max_active_mc = config$max_active_mc,
        active_mc_headroom = config$active_mc_headroom,
        cv_folds = config$cv_folds,
        cv_nlambda = config$cv_nlambda,
        cv_workers_per_fit = config$cv_workers_per_fit,
        backend = config$scheduler_backend,
        legacy_fit_backend = config$backend,
        use_cache = use_cache,
        resumed = FALSE,
        checkpoint_file = checkpoint_file,
        stringsAsFactors = FALSE
      )

      raw_fit_timings <- if (!is.null(assembled)) {
        assembled$fit_timings
      } else {
        valid <- vapply(
          fit_outputs,
          function(x) is.list(x) && is.data.frame(x$timing),
          logical(1)
        )
        bind_rows(lapply(fit_outputs[valid], `[[`, "timing"))
      }

      run_results <- data.frame()
      run_fit_timings <- data.frame()
      if (!is.null(assembled)) {
        run_results <- assembled$results %>%
          mutate(
            run_id = run_id,
            sample_design = sample_design,
            implementation = REFERENCE_STUDY_IMPLEMENTATION,
            cpu_plan = cpu_plan_label,
            scheduler = config$scheduler,
            scheduler_backend = config$scheduler_backend,
            run_index = ii,
            total_runs = nrow(sizes),
            mc_rep = mc_rep_ii,
            mc_reps = mc_reps,
            status = status,
            elapsed_sec = elapsed,
            elapsed_min = elapsed / 60,
            total_cores = config$total_cores,
            reserve_cores = config$reserve_cores,
            usable_cores = config$usable_cores,
            fit_workers = config$fit_workers,
            fit_workers_per_mc = config$fit_workers_per_mc,
            max_concurrent_hal_fits = config$max_concurrent_hal_fits,
            max_active_mc = config$max_active_mc,
            active_mc_headroom = config$active_mc_headroom,
            cv_folds = config$cv_folds,
            cv_nlambda = config$cv_nlambda,
            cv_workers_per_fit = config$cv_workers_per_fit,
            backend = config$scheduler_backend,
            legacy_fit_backend = config$backend,
            use_cache = use_cache,
            resumed = FALSE,
            checkpoint_file = checkpoint_file
          )
      }

      if (nrow(raw_fit_timings) > 0L) {
        run_fit_timings <- raw_fit_timings %>%
          mutate(
            run_id = run_id,
            sample_design = sample_design,
            implementation = REFERENCE_STUDY_IMPLEMENTATION,
            cpu_plan = cpu_plan_label,
            scheduler = config$scheduler,
            scheduler_backend = config$scheduler_backend,
            run_index = ii,
            total_runs = nrow(sizes),
            mc_rep = mc_rep_ii,
            mc_reps = mc_reps,
            seed = seed,
            n_trial = n_trial_ii,
            n_aux = n_aux_ii,
            run_elapsed_sec = elapsed,
            total_cores = config$total_cores,
            reserve_cores = config$reserve_cores,
            usable_cores = config$usable_cores,
            fit_workers = config$fit_workers,
            fit_workers_per_mc = config$fit_workers_per_mc,
            max_concurrent_hal_fits = config$max_concurrent_hal_fits,
            max_active_mc = config$max_active_mc,
            active_mc_headroom = config$active_mc_headroom,
            cv_workers_per_fit = config$cv_workers_per_fit,
            backend = config$scheduler_backend,
            legacy_fit_backend = config$backend,
            use_cache = use_cache,
            resumed = FALSE,
            checkpoint_file = checkpoint_file
          )
      }

      timings[[ii]] <<- timing_row
      fit_timings_all[[ii]] <<- run_fit_timings
      results[[ii]] <<- run_results

      checkpoint <- list(
        schema = RUN_CHECKPOINT_SCHEMA,
        run_signature = run_signature,
        attempt_id = attempt_id,
        status = status,
        timing = timing_row,
        fit_timings = run_fit_timings,
        results = run_results,
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
      )
      atomic_save_rds(
        checkpoint,
        checkpoint_file,
        compress = FALSE,
        overwrite = FALSE
      )
      shard_paths <- write_run_output_shards(
        checkpoint,
        checkpoint_dir,
        ii,
        run_signature,
        attempt_id = attempt_id
      )

      cat(sprintf(
        "[%s] finish ii=%d/%d status=%s elapsed=%.2f sec checkpoint=%s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        ii,
        nrow(sizes),
        status,
        elapsed,
        basename(checkpoint_file)
      ))
      cat("  timing shard = ", shard_paths$timing, "\n", sep = "")
      if (identical(status, "error")) {
        warning("Run ", ii, " failed: ", err_msg, call. = FALSE)
      }
      flush.console()

      list(
        run_index = ii,
        status = status,
        checkpoint_file = checkpoint_file,
        shard_paths = shard_paths
      )
    }

    if (length(pending_indices) > 0L) {
      scheduler_result <- run_global_hal_scheduler(
        run_ids = pending_indices,
        prepare_run = prepare_global_run,
        finalize_run = finalize_global_run,
        project_dir = project_dir,
        global_workers = config$max_concurrent_hal_fits,
        per_run_workers = config$fit_workers_per_mc,
        max_active_runs = config$max_active_mc,
        use_cache = use_cache,
        cv_workers_per_fit = config$cv_workers_per_fit,
        verbose = TRUE
      )
      scheduler_trace <- scheduler_result$trace
      scheduler_stats <- scheduler_result[c(
        "max_observed_global_fits",
        "max_observed_active_runs",
        "max_observed_per_run_fits"
      )]
      finalizer_failed <- vapply(
        scheduler_result$finalized,
        function(x) is.list(x) && isTRUE(x$finalizer_failed),
        logical(1)
      )
      if (any(finalizer_failed)) {
        if (nrow(scheduler_trace) > 0L) {
          atomic_write_csv(
            scheduler_trace,
            scheduler_trace_file,
            overwrite = TRUE
          )
        }
        failed_messages <- vapply(
          scheduler_result$finalized[finalizer_failed],
          function(x) paste0("run ", x$run_id, ": ", x$error),
          character(1)
        )
        stop(
          "One or more run finalizers failed; aggregate output was not ",
          "reported as successful. ",
          paste(failed_messages, collapse = "; "),
          call. = FALSE
        )
      }
    } else if (length(resume_trace_rows) > 0L) {
      scheduler_trace <- bind_rows(resume_trace_rows)
    }
  } else {

  fit_cluster <- NULL
  cv_cluster <- NULL
  worker_pools_started <- FALSE
  on.exit({
    if (!is.null(cv_cluster)) {
      try(parallel::stopCluster(cv_cluster), silent = TRUE)
    }
    if (!is.null(fit_cluster)) {
      try(parallel::stopCluster(fit_cluster), silent = TRUE)
    }
  }, add = TRUE)

  start_worker_pools <- function() {
    if (isTRUE(worker_pools_started)) {
      return(invisible(NULL))
    }

    if (config$fit_workers > 1L &&
        identical(config$backend, "psock")) {
      fit_cluster <<- make_hal_fit_psock_cluster(config$fit_workers)
      cat(
        "Created reusable fit PSOCK pool with ",
        config$fit_workers,
        " workers.\n",
        sep = ""
      )
    } else if (config$fit_workers == 1L &&
               config$cv_workers_per_fit > 1L) {
      cv_cluster <<- make_hal_cv_psock_cluster(
        config$cv_workers_per_fit
      )
      cat(
        "Created reusable CV PSOCK pool with ",
        config$cv_workers_per_fit,
        " workers.\n",
        sep = ""
      )
    }

    worker_pools_started <<- TRUE
    flush.console()
    invisible(NULL)
  }

  for (ii in seq_len(nrow(sizes))) {
    seed <- sizes$seed[[ii]]
    mc_rep_ii <- sizes$mc_rep[[ii]]
    n_trial_ii <- sizes$n_trial[[ii]]
    n_aux_ii <- sizes$n_aux[[ii]]
    run_signature <- make_run_signature(
      project_dir = project_dir,
      seed = seed,
      m = m,
      n_trial = n_trial_ii,
      n_aux = n_aux_ii,
      p_resp = p_resp,
      p_cens = p_cens,
      config = config,
      code_fingerprint = code_fingerprint
    )
    checkpoint_match <- if (isTRUE(resume)) {
      find_successful_run_checkpoint(
        checkpoint_dir,
        ii,
        run_signature
      )
    } else {
      list(checkpoint = NULL, path = NULL)
    }
    checkpoint <- checkpoint_match$checkpoint

    if (!is.null(checkpoint)) {
      checkpoint_file <- checkpoint_match$path
      write_run_output_shards(
        checkpoint,
        checkpoint_dir,
        ii,
        run_signature,
        attempt_id = checkpoint$attempt_id
      )
      checkpoint$timing <- refresh_resumed_metadata(
        checkpoint$timing,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )
      checkpoint$fit_timings <- refresh_resumed_metadata(
        checkpoint$fit_timings,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )
      checkpoint$results <- refresh_resumed_metadata(
        checkpoint$results,
        run_id,
        sample_design,
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        cpu_plan_label,
        config,
        use_cache,
        checkpoint_file
      )

      timings[[ii]] <- checkpoint$timing
      fit_timings_all[[ii]] <- checkpoint$fit_timings
      results[[ii]] <- checkpoint$results
      cat(sprintf(
        "[%s] resume ii=%d/%d mc=%d/%d seed=%d from %s\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        ii,
        nrow(sizes),
        mc_rep_ii,
        mc_reps,
        seed,
        basename(checkpoint_file)
      ))
      flush.console()
      next
    }

    cat(sprintf(
      "[%s] start ii=%d/%d mc=%d/%d seed=%d n_trial=%d n_aux=%d\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ii,
      nrow(sizes),
      mc_rep_ii,
      mc_reps,
      seed,
      n_trial_ii,
      n_aux_ii
    ))
    flush.console()

    gc()
    t0 <- proc.time()
    status <- "ok"
    err_msg <- NA_character_

    result <- tryCatch(
      {
        start_worker_pools()
        run_reference_replicate(
          m = m,
          n_trial = n_trial_ii,
          n_aux = n_aux_ii,
          p_resp = p_resp,
          p_cens = p_cens,
          seed = seed,
          config = config,
          cache_dir = cache_dir,
          use_cache = use_cache,
          fit_cluster = fit_cluster,
          cv_cluster = cv_cluster
        )
      },
      error = function(e) {
        status <<- "error"
        err_msg <<- conditionMessage(e)
        NULL
      }
    )

    elapsed <- unname((proc.time() - t0)[["elapsed"]])
    attempt_id <- make_run_attempt_id()
    checkpoint_file <- make_run_checkpoint_file(
      checkpoint_dir,
      ii,
      run_signature,
      attempt_id = attempt_id
    )
    timing_row <- data.frame(
      run_id = run_id,
      sample_design = sample_design,
      implementation = REFERENCE_STUDY_IMPLEMENTATION,
      cpu_plan = cpu_plan_label,
      run_index = ii,
      total_runs = nrow(sizes),
      mc_rep = mc_rep_ii,
      mc_reps = mc_reps,
      seed = seed,
      n_trial = n_trial_ii,
      n_aux = n_aux_ii,
      status = status,
      elapsed_sec = elapsed,
      elapsed_min = elapsed / 60,
      result_rows = if (is.null(result)) NA_integer_ else nrow(result$results),
      error = err_msg,
      total_cores = config$total_cores,
      reserve_cores = config$reserve_cores,
      usable_cores = config$usable_cores,
      fit_workers = config$fit_workers,
      cv_folds = config$cv_folds,
      cv_nlambda = config$cv_nlambda,
      cv_workers_per_fit = config$cv_workers_per_fit,
      backend = config$backend,
      use_cache = use_cache,
      resumed = FALSE,
      checkpoint_file = checkpoint_file,
      stringsAsFactors = FALSE
    )
    timings[[ii]] <- timing_row

    run_results <- data.frame()
    run_fit_timings <- data.frame()
    if (!is.null(result)) {
      run_results <- result$results %>%
        mutate(
          run_id = run_id,
          sample_design = sample_design,
          implementation = REFERENCE_STUDY_IMPLEMENTATION,
          cpu_plan = cpu_plan_label,
          run_index = ii,
          total_runs = nrow(sizes),
          mc_rep = mc_rep_ii,
          mc_reps = mc_reps,
          status = status,
          elapsed_sec = elapsed,
          elapsed_min = elapsed / 60,
          total_cores = config$total_cores,
          reserve_cores = config$reserve_cores,
          usable_cores = config$usable_cores,
          fit_workers = config$fit_workers,
          cv_folds = config$cv_folds,
          cv_nlambda = config$cv_nlambda,
          cv_workers_per_fit = config$cv_workers_per_fit,
          backend = config$backend,
          use_cache = use_cache,
          resumed = FALSE,
          checkpoint_file = checkpoint_file
        )

      run_fit_timings <- result$fit_timings %>%
        mutate(
          run_id = run_id,
          sample_design = sample_design,
          implementation = REFERENCE_STUDY_IMPLEMENTATION,
          cpu_plan = cpu_plan_label,
          run_index = ii,
          total_runs = nrow(sizes),
          mc_rep = mc_rep_ii,
          mc_reps = mc_reps,
          seed = seed,
          n_trial = n_trial_ii,
          n_aux = n_aux_ii,
          run_elapsed_sec = elapsed,
          total_cores = config$total_cores,
          reserve_cores = config$reserve_cores,
          usable_cores = config$usable_cores,
          fit_workers = config$fit_workers,
          cv_workers_per_fit = config$cv_workers_per_fit,
          backend = config$backend,
          use_cache = use_cache,
          resumed = FALSE,
          checkpoint_file = checkpoint_file
        )

      results[[ii]] <- run_results
      fit_timings_all[[ii]] <- run_fit_timings
    }

    checkpoint <- list(
      schema = RUN_CHECKPOINT_SCHEMA,
      run_signature = run_signature,
      attempt_id = attempt_id,
      status = status,
      timing = timing_row,
      fit_timings = run_fit_timings,
      results = run_results,
      created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
    atomic_save_rds(
      checkpoint,
      checkpoint_file,
      compress = FALSE,
      overwrite = FALSE
    )
    shard_paths <- write_run_output_shards(
      checkpoint,
      checkpoint_dir,
      ii,
      run_signature,
      attempt_id = attempt_id
    )

    if (identical(status, "error")) {
      if (!is.null(cv_cluster)) {
        try(parallel::stopCluster(cv_cluster), silent = TRUE)
        cv_cluster <- NULL
      }
      if (!is.null(fit_cluster)) {
        try(parallel::stopCluster(fit_cluster), silent = TRUE)
        fit_cluster <- NULL
      }
      worker_pools_started <- FALSE
      warning(
        "Run ", ii, " failed; reusable worker pools were reset before ",
        "the next run. Error: ", err_msg,
        call. = FALSE
      )
    }

    cat(sprintf(
      "[%s] finish ii=%d/%d status=%s elapsed=%.2f sec checkpoint=%s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ii,
      nrow(sizes),
      status,
      elapsed,
      basename(checkpoint_file)
    ))
    cat("  timing shard = ", shard_paths$timing, "\n", sep = "")
    flush.console()
  }
  }

  timing_output <- bind_rows(timings)
  fit_timing_output <- bind_rows(fit_timings_all)
  result_output <- bind_rows(results)

  atomic_write_csv(timing_output, timing_file, overwrite = TRUE)
  if (nrow(scheduler_trace) > 0L) {
    atomic_write_csv(
      scheduler_trace,
      scheduler_trace_file,
      overwrite = TRUE
    )
  } else if (file.exists(scheduler_trace_file)) {
    quarantine_invalid_file(scheduler_trace_file, "no-current-scheduler-events")
  }
  if (nrow(fit_timing_output) > 0L) {
    atomic_write_csv(
      fit_timing_output,
      fit_timing_file,
      overwrite = TRUE
    )
  } else if (file.exists(fit_timing_file)) {
    quarantine_invalid_file(fit_timing_file, "no-current-fit-results")
  }
  if (nrow(result_output) > 0L) {
    atomic_write_csv(result_output, result_file, overwrite = TRUE)
  } else {
    if (file.exists(result_file)) {
      quarantine_invalid_file(result_file, "no-current-results")
    }
  }

  invocation_end_time <- Sys.time()
  invocation_wall_elapsed_sec <- as.numeric(difftime(
    invocation_end_time,
    invocation_start_time,
    units = "secs"
  ))
  complete_run_rows <- nrow(timing_output) == nrow(sizes) &&
    setequal(as.integer(timing_output$run_index), seq_len(nrow(sizes)))
  invocation_status <- if (isTRUE(complete_run_rows) &&
                           all(as.character(timing_output$status) == "ok")) {
    "ok"
  } else {
    "error"
  }
  run_summary <- data.frame(
    run_id = run_id,
    sample_design = sample_design,
    implementation = REFERENCE_STUDY_IMPLEMENTATION,
    status = invocation_status,
    sample_size_combinations = length(unique(paste(
      sizes$n_trial,
      sizes$n_aux,
      sep = "/"
    ))),
    mc_reps = mc_reps,
    total_runs = nrow(sizes),
    successful_runs = sum(
      as.character(timing_output$status) == "ok",
      na.rm = TRUE
    ),
    resumed_runs = sum(timing_output$resumed %in% TRUE, na.rm = TRUE),
    estimator_rows = nrow(result_output),
    fit_timing_rows = nrow(fit_timing_output),
    invocation_start_time = format(
      invocation_start_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    invocation_end_time = format(
      invocation_end_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    invocation_wall_elapsed_sec = invocation_wall_elapsed_sec,
    invocation_wall_elapsed_min = invocation_wall_elapsed_sec / 60,
    invocation_wall_elapsed_hour = invocation_wall_elapsed_sec / 3600,
    recorded_run_elapsed_sec = sum(
      timing_output$elapsed_sec,
      na.rm = TRUE
    ),
    total_cores = config$total_cores,
    reserve_cores = config$reserve_cores,
    fit_workers = config$fit_workers,
    fit_workers_per_mc = config$fit_workers_per_mc,
    cv_workers_per_fit = config$cv_workers_per_fit,
    max_concurrent_hal_fits = config$max_concurrent_hal_fits,
    max_active_mc = config$max_active_mc,
    active_mc_headroom = config$active_mc_headroom,
    max_observed_global_fits = scheduler_stats$max_observed_global_fits,
    max_observed_active_runs = scheduler_stats$max_observed_active_runs,
    max_observed_per_run_fits = scheduler_stats$max_observed_per_run_fits,
    cv_folds = config$cv_folds,
    cv_nlambda = config$cv_nlambda,
    scheduler = config$scheduler,
    scheduler_backend = config$scheduler_backend,
    backend = config$scheduler_backend,
    legacy_fit_backend = config$backend,
    use_cache = isTRUE(use_cache),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(run_summary, summary_file, overwrite = TRUE)

  cat("Finished the PopART manuscript-reference Monte Carlo study\n")
  cat("Wrote parallel plan to: ", plan_file, "\n", sep = "")
  cat("Wrote timing results to: ", timing_file, "\n", sep = "")
  if (nrow(scheduler_trace) > 0L) {
    cat("Wrote scheduler trace to: ", scheduler_trace_file, "\n", sep = "")
  }
  if (nrow(fit_timing_output) > 0L) {
    cat("Wrote fit timing results to: ", fit_timing_file, "\n", sep = "")
  }
  if (nrow(result_output) > 0L) {
    cat("Wrote estimator results to: ", result_file, "\n", sep = "")
  }
  cat("Wrote invocation summary to: ", summary_file, "\n", sep = "")
  cat("Checkpoint and CSV shards are in: ", checkpoint_dir, "\n", sep = "")

  invisible(list(
    config = config,
    sample_grid = sizes,
    files = output_files,
    timings = timing_output,
    fit_timings = fit_timing_output,
    results = result_output,
    run_summary = run_summary,
    scheduler_trace = scheduler_trace,
    scheduler_stats = scheduler_stats
  ))
}
