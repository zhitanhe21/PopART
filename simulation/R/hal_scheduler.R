###############################################################################
# Cross-replicate global HAL-fit queue
###############################################################################


# Mutable state owned by one R process. In an installed package this
# environment lives inside the package namespace; in source mode it lives in
# the environment into which the helpers were sourced. Keeping the private CV
# pool here avoids relying on objects injected into `.GlobalEnv`.
.aipw_global_worker_state <- new.env(parent = emptyenv())


#' Initialize process-local resources for one persistent HAL slot.
#'
#' @param cv_workers_per_fit Positive integer number of CV workers owned by the
#'   complete-HAL slot.
#'
#' @return A list containing the slot PID and effective CV worker count.
initialize_global_hal_worker <- function(cv_workers_per_fit = 1L) {
  cv_workers_per_fit <- suppressWarnings(as.integer(cv_workers_per_fit))
  if (is.na(cv_workers_per_fit) || cv_workers_per_fit < 1L) {
    stop("cv_workers_per_fit must be a positive integer.", call. = FALSE)
  }

  # Reinitialization must never leak an earlier private pool.
  shutdown_global_hal_worker()
  set_single_threaded_math()
  foreach::registerDoSEQ()

  if (cv_workers_per_fit > 1L) {
    .aipw_global_worker_state$cv_cluster <-
      make_hal_cv_psock_cluster(cv_workers_per_fit)
  } else {
    .aipw_global_worker_state$cv_cluster <- NULL
  }
  .aipw_global_worker_state$cv_workers_per_fit <- cv_workers_per_fit

  list(pid = Sys.getpid(), cv_workers = cv_workers_per_fit)
}


#' Release process-local resources owned by one persistent HAL slot.
#'
#' @return Invisibly returns `TRUE`.
shutdown_global_hal_worker <- function() {
  if (exists(
    "cv_cluster",
    envir = .aipw_global_worker_state,
    inherits = FALSE
  )) {
    cv_cluster <- get(
      "cv_cluster",
      envir = .aipw_global_worker_state,
      inherits = FALSE
    )
    if (!is.null(cv_cluster)) {
      try(parallel::stopCluster(cv_cluster), silent = TRUE)
    }
    rm("cv_cluster", envir = .aipw_global_worker_state)
  }
  if (exists(
    "cv_workers_per_fit",
    envir = .aipw_global_worker_state,
    inherits = FALSE
  )) {
    rm("cv_workers_per_fit", envir = .aipw_global_worker_state)
  }
  if (requireNamespace("foreach", quietly = TRUE)) {
    foreach::registerDoSEQ()
  }
  invisible(TRUE)
}


#' Execute one task inside a persistent callr HAL slot.
#'
#' The wrapper always returns a serializable success/error envelope. Ordinary
#' HAL or glmnet errors therefore fail one fit without terminating the slot
#' process or unrelated runs.
#' @noRd
global_hal_task_envelope <- function(task_fun, spec, use_cache, task_id) {
  started_at <- Sys.time()
  tryCatch(
    list(
      ok = TRUE,
      value = task_fun(spec, use_cache),
      error = NA_character_,
      task_id = task_id,
      pid = Sys.getpid(),
      started_at = format(started_at, "%Y-%m-%d %H:%M:%OS6"),
      finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6")
    ),
    error = function(e) {
      list(
        ok = FALSE,
        value = NULL,
        error = conditionMessage(e),
        error_class = class(e),
        task_id = task_id,
        pid = Sys.getpid(),
        started_at = format(started_at, "%Y-%m-%d %H:%M:%OS6"),
        finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6")
      )
    }
  )
}


#' Default complete-HAL task used by a persistent global slot.
#' @noRd
global_hal_default_task <- function(spec, use_cache) {
  function_env <- environment(global_hal_default_task)
  worker_state <- get0(
    ".aipw_global_worker_state",
    envir = function_env,
    inherits = TRUE,
    ifnotfound = NULL
  )
  cv_cluster <- if (is.environment(worker_state)) {
    get0(
      "cv_cluster",
      envir = worker_state,
      inherits = FALSE,
      ifnotfound = NULL
    )
  } else {
    # Compatibility with slots initialized by an older sourced implementation.
    get0(
      ".AIPW_GLOBAL_CV_CLUSTER",
      envir = .GlobalEnv,
      inherits = FALSE,
      ifnotfound = NULL
    )
  }
  run_job <- get0(
    "run_one_hal_job",
    envir = function_env,
    inherits = TRUE,
    ifnotfound = NULL
  )
  if (!is.function(run_job) && isNamespaceLoaded("popart")) {
    run_job <- get("run_one_hal_job", envir = asNamespace("popart"))
  }
  if (!is.function(run_job)) {
    stop("run_one_hal_job() is unavailable in this HAL slot.", call. = FALSE)
  }
  run_job(
    spec,
    use_cache = use_cache,
    cv_cluster = cv_cluster
  )
}


#' Create persistent callr sessions representing global HAL slots.
#'
#' With `cv_workers_per_fit = 1`, each session directly runs one complete HAL
#' fit. With larger values, each session owns a private reusable CV PSOCK pool;
#' the socket cluster never crosses a process boundary.
#' @noRd
make_global_hal_callr_slots <- function(worker_count, project_dir = NULL,
                                        cv_workers_per_fit = 1L,
                                        initialize_hal = TRUE) {
  if (!requireNamespace("callr", quietly = TRUE)) {
    stop(
      "The global scheduler requires the callr package. Install callr and retry.",
      call. = FALSE
    )
  }

  worker_count <- suppressWarnings(as.integer(worker_count))
  cv_workers_per_fit <- suppressWarnings(as.integer(cv_workers_per_fit))
  if (is.na(worker_count) || worker_count < 1L) {
    stop("worker_count must be a positive integer.", call. = FALSE)
  }
  if (is.na(cv_workers_per_fit) || cv_workers_per_fit < 1L) {
    stop("cv_workers_per_fit must be a positive integer.", call. = FALSE)
  }

  source_mode <- FALSE
  package_name <- "popart"
  if (isTRUE(initialize_hal)) {
    source_file <- if (!is.null(project_dir) && length(project_dir) == 1L &&
                       !is.na(project_dir) && nzchar(project_dir)) {
      file.path(project_dir, "simulation", "R", "load_simulation.R")
    } else {
      NA_character_
    }
    source_mode <- !is.na(source_file) && file.exists(source_file)
    if (isTRUE(source_mode)) {
      project_dir <- normalizePath(project_dir)
    } else if (!requireNamespace(package_name, quietly = TRUE)) {
      stop(
        paste0(
          "HAL slot initialization requires either a source project_dir ",
          "containing simulation/R/load_simulation.R or an installed popart package."
        ),
        call. = FALSE
      )
    }
  }

  slots <- vector("list", worker_count)
  initialized <- 0L
  creation_complete <- FALSE
  on.exit({
    if (!isTRUE(creation_complete)) {
      close_global_hal_callr_slots(slots[seq_len(max(0L, initialized))])
    }
  }, add = TRUE)

  for (ii in seq_len(worker_count)) {
    options <- callr::r_session_options(
      supervise = TRUE,
      system_profile = FALSE,
      user_profile = FALSE,
      libpath = .libPaths(),
      env = c(
        TERM = "dumb",
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
    )
    slot <- callr::r_session$new(
      options = options,
      wait = TRUE,
      wait_timeout = 10000
    )
    # Register the session for cleanup before project initialization. Package
    # loading or private CV-pool creation can fail after the process starts.
    slots[[ii]] <- slot
    initialized <- ii

    if (isTRUE(initialize_hal)) {
      slot$run(
        function(project_dir, source_mode, package_name,
                 cv_workers_per_fit) {
          Sys.setenv(
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            NUMEXPR_NUM_THREADS = "1"
          )

          if (isTRUE(source_mode)) {
            setwd(project_dir)
            source(file.path(
              project_dir, "simulation", "R", "load_simulation.R"
            ))
            source_reference_study(project_dir, envir = .GlobalEnv)
            load_reference_study_packages()
            initialize_worker <- get(
              "initialize_global_hal_worker",
              envir = .GlobalEnv,
              inherits = TRUE
            )
          } else {
            namespace <- loadNamespace(package_name)
            initialize_worker <- get(
              "initialize_global_hal_worker",
              envir = namespace,
              inherits = FALSE
            )
          }
          initialize_worker(cv_workers_per_fit)
        },
        args = list(
          project_dir = project_dir,
          source_mode = source_mode,
          package_name = package_name,
          cv_workers_per_fit = cv_workers_per_fit
        )
      )
    }
  }

  creation_complete <- TRUE
  slots
}


#' Close persistent global HAL slots and their private CV pools.
#' @noRd
close_global_hal_callr_slots <- function(slots, grace = 1000L) {
  if (length(slots) == 0L) return(invisible(NULL))
  for (slot in slots) {
    if (is.null(slot)) next
    state <- tryCatch(slot$get_state(), error = function(e) "finished")
    if (identical(state, "idle")) {
      try(
        slot$run(function() {
          shutdown_worker <- get0(
            "shutdown_global_hal_worker",
            envir = .GlobalEnv,
            inherits = TRUE,
            ifnotfound = NULL
          )
          if (!is.function(shutdown_worker) &&
              isNamespaceLoaded("popart")) {
            shutdown_worker <- get(
              "shutdown_global_hal_worker",
              envir = asNamespace("popart"),
              inherits = FALSE
            )
          }
          if (is.function(shutdown_worker)) {
            shutdown_worker()
          } else {
            # Compatibility cleanup for older sourced slot initialization.
            cv_cluster <- get0(
              ".AIPW_GLOBAL_CV_CLUSTER",
              envir = .GlobalEnv,
              inherits = FALSE,
              ifnotfound = NULL
            )
            if (!is.null(cv_cluster)) {
              try(parallel::stopCluster(cv_cluster), silent = TRUE)
              rm(".AIPW_GLOBAL_CV_CLUSTER", envir = .GlobalEnv)
            }
            if (requireNamespace("foreach", quietly = TRUE)) {
              foreach::registerDoSEQ()
            }
          }
          TRUE
        }),
        silent = TRUE
      )
    }
    try(slot$close(grace = grace), silent = TRUE)
  }
  invisible(NULL)
}


read_global_hal_slot_event <- function(slot) {
  repeat {
    event <- tryCatch(slot$read(), error = function(e) e)
    if (inherits(event, "error")) {
      return(list(ok = FALSE, error = conditionMessage(event), crashed = TRUE))
    }
    if (is.null(event)) return(NULL)

    code <- suppressWarnings(as.integer(event$code))
    if (identical(code, 200L)) {
      if (!is.null(event$error)) {
        return(list(
          ok = FALSE,
          error = conditionMessage(event$error),
          crashed = FALSE
        ))
      }
      return(list(ok = TRUE, envelope = event$result, crashed = FALSE))
    }
    if (code %in% c(500L, 501L, 502L)) {
      return(list(
        ok = FALSE,
        error = if (!is.null(event$message)) {
          as.character(event$message)
        } else {
          paste("callr slot terminated with event code", code)
        },
        crashed = TRUE
      ))
    }
    # Message/progress events do not complete the task; keep draining.
  }
}


#' Run complete HAL fits from multiple MC runs through one shared queue.
#'
#' @param run_ids Unique run identifiers in admission order.
#' @param prepare_run Function taking one run ID and returning a list with named
#'   `jobs`; every other component is retained as the finalization context.
#' @param finalize_run Function called as
#'   `finalize_run(run_id, context, fit_outputs, error, metrics)`.
#' @param project_dir Optional repository source directory used to initialize
#'   HAL slots. When omitted outside source mode, workers load the installed
#'   the `popart` namespace.
#' @param global_workers Maximum concurrently running complete HAL fits.
#' @param cv_workers_per_fit Positive integer number of CV workers owned by
#'   each complete-HAL slot.
#' @param per_run_workers Maximum concurrently running fits from one run.
#' @param max_active_runs Maximum prepared, unfinished run contexts.
#' @param use_cache Whether workers may use the fit cache.
#' @param task_fun Function executed for each job inside a slot.
#' @param slots Optional pre-created callr sessions, mainly for tests.
#' @param initialize_hal Whether owned slots should load the HAL project.
#' @param poll_interval_ms Milliseconds between completion polls.
#' @param verbose Whether to print admission/completion progress.
#'
#' @return A list with finalizer return values, an event trace, and observed
#'   concurrency maxima.
#' @export
run_global_hal_scheduler <- function(
    run_ids,
    prepare_run,
    finalize_run,
    project_dir = NULL,
    global_workers = 1L,
    per_run_workers = 2L,
    max_active_runs = 1L,
    use_cache = TRUE,
    task_fun = global_hal_default_task,
    slots = NULL,
    initialize_hal = TRUE,
    cv_workers_per_fit = 1L,
    poll_interval_ms = 50L,
    verbose = TRUE) {

  if (!is.function(prepare_run) || !is.function(finalize_run) ||
      !is.function(task_fun)) {
    stop("prepare_run, finalize_run, and task_fun must be functions.", call. = FALSE)
  }
  if (length(run_ids) < 1L || anyDuplicated(as.character(run_ids))) {
    stop("run_ids must be nonempty and unique.", call. = FALSE)
  }

  global_workers <- suppressWarnings(as.integer(global_workers))
  per_run_workers <- suppressWarnings(as.integer(per_run_workers))
  max_active_runs <- suppressWarnings(as.integer(max_active_runs))
  poll_interval_ms <- suppressWarnings(as.integer(poll_interval_ms))
  if (is.na(global_workers) || global_workers < 1L ||
      is.na(per_run_workers) || per_run_workers < 1L ||
      is.na(max_active_runs) || max_active_runs < 1L) {
    stop("Scheduler limits must be positive integers.", call. = FALSE)
  }
  poll_interval_ms <- min(1000L, max(10L, poll_interval_ms))

  global_workers <- min(global_workers, length(run_ids) * per_run_workers)
  owns_slots <- is.null(slots)
  if (isTRUE(owns_slots)) {
    slots <- make_global_hal_callr_slots(
      global_workers,
      project_dir = project_dir,
      cv_workers_per_fit = cv_workers_per_fit,
      initialize_hal = initialize_hal
    )
  } else if (length(slots) < global_workers) {
    stop("The supplied slot list is smaller than global_workers.", call. = FALSE)
  } else {
    slots <- slots[seq_len(global_workers)]
  }
  if (isTRUE(owns_slots)) {
    on.exit(close_global_hal_callr_slots(slots), add = TRUE)
  }

  scheduler_started <- Sys.time()
  states <- new.env(hash = TRUE, parent = emptyenv())
  active_keys <- character()
  next_run <- 1L
  round_robin_cursor <- 0L
  slot_busy <- rep(FALSE, global_workers)
  slot_meta <- vector("list", global_workers)
  slot_alive <- rep(TRUE, global_workers)
  task_counter <- 0L
  finalized <- list()
  trace_rows <- list()
  trace_counter <- 0L
  max_observed_global <- 0L
  max_observed_active <- 0L
  max_observed_per_run <- 0L

  trace_event <- function(event, run_id = NA, fit_label = NA_character_,
                          node = NA_integer_, task_id = NA_character_,
                          message = NA_character_) {
    trace_counter <<- trace_counter + 1L
    running_global <- sum(slot_busy)
    per_run_running <- if (length(active_keys) == 0L) {
      0L
    } else {
      max(vapply(
        active_keys,
        function(key) get(key, envir = states, inherits = FALSE)$running,
        integer(1)
      ))
    }
    max_observed_global <<- max(max_observed_global, running_global)
    max_observed_active <<- max(max_observed_active, length(active_keys))
    max_observed_per_run <<- max(max_observed_per_run, per_run_running)
    trace_rows[[trace_counter]] <<- data.frame(
      event_index = trace_counter,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%OS6"),
      elapsed_sec = as.numeric(difftime(Sys.time(), scheduler_started, units = "secs")),
      event = event,
      run_id = as.character(run_id),
      fit_label = fit_label,
      node = as.integer(node),
      task_id = task_id,
      active_runs = length(active_keys),
      running_fits = running_global,
      free_slots = sum(!slot_busy & slot_alive),
      max_running_this_run = per_run_running,
      message = message,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  make_metrics <- function(state, finished_at = Sys.time()) {
    list(
      admitted_at = state$admitted_at,
      prepared_at = state$prepared_at,
      first_dispatched_at = state$first_dispatched_at,
      last_completed_at = state$last_completed_at,
      prepare_elapsed_sec = as.numeric(difftime(
        state$prepared_at,
        state$admitted_at,
        units = "secs"
      )),
      queue_wait_sec = if (is.na(state$first_dispatched_at)) {
        NA_real_
      } else {
        as.numeric(difftime(
          state$first_dispatched_at,
          state$prepared_at,
          units = "secs"
        ))
      },
      fit_span_sec = if (is.na(state$first_dispatched_at) ||
                         is.na(state$last_completed_at)) {
        NA_real_
      } else {
        as.numeric(difftime(
          state$last_completed_at,
          state$first_dispatched_at,
          units = "secs"
        ))
      },
      run_wall_elapsed_sec = as.numeric(difftime(
        finished_at,
        state$admitted_at,
        units = "secs"
      ))
    )
  }

  safely_finalize <- function(run_id, context, fit_outputs, error, metrics) {
    tryCatch(
      {
        value <- finalize_run(run_id, context, fit_outputs, error, metrics)
        if (!is.list(value)) {
          stop("finalize_run must return a non-NULL list.", call. = FALSE)
        }
        value
      },
      error = function(e) {
        list(
          run_id = run_id,
          status = "error",
          finalizer_failed = TRUE,
          error = paste("finalize_run failed:", conditionMessage(e))
        )
      }
    )
  }

  admit_runs <- function() {
    while (length(active_keys) < max_active_runs && next_run <= length(run_ids)) {
      run_id <- run_ids[[next_run]]
      next_run <<- next_run + 1L
      key <- as.character(run_id)
      admitted_at <- Sys.time()
      prepared <- tryCatch(prepare_run(run_id), error = function(e) e)

      if (inherits(prepared, "error")) {
        error_message <- paste("run preparation failed:", conditionMessage(prepared))
        metrics <- list(
          admitted_at = admitted_at,
          prepared_at = Sys.time(),
          first_dispatched_at = as.POSIXct(NA),
          last_completed_at = as.POSIXct(NA),
          prepare_elapsed_sec = as.numeric(difftime(
            Sys.time(), admitted_at, units = "secs"
          )),
          queue_wait_sec = NA_real_,
          fit_span_sec = NA_real_,
          run_wall_elapsed_sec = as.numeric(difftime(
            Sys.time(), admitted_at, units = "secs"
          ))
        )
        trace_event("prepare_error", run_id = run_id, message = error_message)
        finalized[[key]] <<- safely_finalize(
          run_id, NULL, list(), error_message, metrics
        )
        next
      }

      jobs <- prepared$jobs
      if (!is.list(jobs) || length(jobs) < 1L || is.null(names(jobs)) ||
          any(!nzchar(names(jobs))) || anyDuplicated(names(jobs))) {
        error_message <- "prepare_run must return uniquely named jobs."
        metrics <- list(
          admitted_at = admitted_at,
          prepared_at = Sys.time(),
          first_dispatched_at = as.POSIXct(NA),
          last_completed_at = as.POSIXct(NA),
          prepare_elapsed_sec = as.numeric(difftime(
            Sys.time(), admitted_at, units = "secs"
          )),
          queue_wait_sec = NA_real_,
          fit_span_sec = NA_real_,
          run_wall_elapsed_sec = as.numeric(difftime(
            Sys.time(), admitted_at, units = "secs"
          ))
        )
        trace_event("prepare_error", run_id = run_id, message = error_message)
        prepared$jobs <- NULL
        finalized[[key]] <<- safely_finalize(
          run_id, prepared, list(), error_message, metrics
        )
        next
      }

      prepared$jobs <- NULL
      now <- Sys.time()
      state <- list(
        run_id = run_id,
        context = prepared,
        jobs = jobs,
        all_labels = names(jobs),
        pending_labels = names(jobs),
        outputs = list(),
        running = 0L,
        error = NULL,
        admitted_at = admitted_at,
        prepared_at = now,
        first_dispatched_at = as.POSIXct(NA),
        last_completed_at = as.POSIXct(NA)
      )
      assign(key, state, envir = states)
      active_keys <<- c(active_keys, key)
      trace_event("admit", run_id = run_id)
      if (isTRUE(verbose)) {
        cat(sprintf(
          "[%s] admit run=%s active=%d/%d\n",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          as.character(run_id),
          length(active_keys),
          max_active_runs
        ))
        flush.console()
      }
    }
    invisible(NULL)
  }

  choose_eligible_run <- function() {
    n_active <- length(active_keys)
    if (n_active == 0L) return(NULL)
    for (step in seq_len(n_active)) {
      round_robin_cursor <<- (round_robin_cursor %% n_active) + 1L
      key <- active_keys[[round_robin_cursor]]
      state <- get(key, envir = states, inherits = FALSE)
      if (is.null(state$error) && state$running < per_run_workers &&
          length(state$pending_labels) > 0L) {
        return(key)
      }
    }
    NULL
  }

  finalize_ready_runs <- function() {
    if (length(active_keys) == 0L) return(invisible(FALSE))
    finalized_any <- FALSE
    for (key in active_keys) {
      state <- get(key, envir = states, inherits = FALSE)
      success_ready <- is.null(state$error) && state$running == 0L &&
        length(state$pending_labels) == 0L &&
        length(state$outputs) == length(state$all_labels)
      error_ready <- !is.null(state$error) && state$running == 0L
      if (!success_ready && !error_ready) next

      fit_outputs <- state$outputs[state$all_labels[
        state$all_labels %in% names(state$outputs)
      ]]
      metrics <- make_metrics(state)
      trace_event(
        if (success_ready) "run_complete" else "run_error",
        run_id = state$run_id,
        message = if (is.null(state$error)) NA_character_ else state$error
      )
      finalized[[key]] <<- safely_finalize(
        state$run_id,
        state$context,
        fit_outputs,
        state$error,
        metrics
      )

      rm(list = key, envir = states)
      active_keys <<- active_keys[active_keys != key]
      round_robin_cursor <<- min(round_robin_cursor, length(active_keys))
      finalized_any <- TRUE
      if (isTRUE(verbose)) {
        cat(sprintf(
          "[%s] finalize run=%s status=%s active=%d\n",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          as.character(state$run_id),
          if (success_ready) "ok" else "error",
          length(active_keys)
        ))
        flush.console()
      }
      rm(state, fit_outputs)
      gc(verbose = FALSE)
    }
    invisible(finalized_any)
  }

  dispatch_jobs <- function() {
    dispatched <- FALSE
    free_nodes <- which(!slot_busy & slot_alive)
    for (node in free_nodes) {
      key <- choose_eligible_run()
      if (is.null(key)) break

      state <- get(key, envir = states, inherits = FALSE)
      label <- state$pending_labels[[1L]]
      spec <- state$jobs[[label]]
      task_counter <<- task_counter + 1L
      task_id <- sprintf("run%s-%s-%06d", key, label, task_counter)
      submitted <- tryCatch(
        {
          slots[[node]]$call(
            global_hal_task_envelope,
            args = list(
              task_fun = task_fun,
              spec = spec,
              use_cache = isTRUE(use_cache),
              task_id = task_id
            ),
            package = FALSE
          )
          TRUE
        },
        error = function(e) e
      )

      if (inherits(submitted, "error")) {
        # A session that rejects an asynchronous call is no longer trusted.
        # Retire the slot but leave the fit pending so another live slot can
        # execute it. If every slot dies, the deadlock guard below fails the
        # remaining runs explicitly instead of repeatedly poisoning them.
        slot_alive[[node]] <<- FALSE
        try(slots[[node]]$close(grace = 100L), silent = TRUE)
        trace_event(
          "submit_error",
          run_id = state$run_id,
          fit_label = label,
          node = node,
          task_id = task_id,
          message = paste(
            "slot retired after task submission failed:",
            conditionMessage(submitted)
          )
        )
        next
      }

      state$pending_labels <- state$pending_labels[-1L]
      state$jobs[[label]] <- NULL
      state$running <- state$running + 1L
      if (is.na(state$first_dispatched_at)) state$first_dispatched_at <- Sys.time()
      assign(key, state, envir = states)
      slot_busy[[node]] <<- TRUE
      slot_meta[[node]] <<- list(
        key = key,
        run_id = state$run_id,
        fit_label = label,
        task_id = task_id,
        submitted_at = Sys.time()
      )
      trace_event(
        "dispatch",
        run_id = state$run_id,
        fit_label = label,
        node = node,
        task_id = task_id
      )
      dispatched <- TRUE
    }
    invisible(dispatched)
  }

  process_slot_result <- function(node, result) {
    meta <- slot_meta[[node]]
    slot_busy[[node]] <<- FALSE
    slot_meta[node] <<- list(NULL)
    if (is.null(meta)) return(invisible(NULL))

    key <- meta$key
    if (!exists(key, envir = states, inherits = FALSE)) {
      trace_event(
        "orphan_result",
        run_id = meta$run_id,
        fit_label = meta$fit_label,
        node = node,
        task_id = meta$task_id
      )
      return(invisible(NULL))
    }
    state <- get(key, envir = states, inherits = FALSE)
    state$running <- max(0L, state$running - 1L)
    state$last_completed_at <- Sys.time()

    envelope_ok <- isTRUE(result$ok) && is.list(result$envelope) &&
      isTRUE(result$envelope$ok)
    if (isTRUE(envelope_ok)) {
      state$outputs[[meta$fit_label]] <- result$envelope$value
      trace_event(
        "complete",
        run_id = meta$run_id,
        fit_label = meta$fit_label,
        node = node,
        task_id = meta$task_id
      )
    } else {
      error_message <- if (isTRUE(result$ok) && is.list(result$envelope)) {
        result$envelope$error
      } else {
        result$error
      }
      if (is.null(error_message) || is.na(error_message) || !nzchar(error_message)) {
        error_message <- "unknown global HAL slot error"
      }
      state$error <- paste0(meta$fit_label, ": ", error_message)
      state$pending_labels <- character()
      trace_event(
        "fit_error",
        run_id = meta$run_id,
        fit_label = meta$fit_label,
        node = node,
        task_id = meta$task_id,
        message = state$error
      )
      if (isTRUE(result$crashed)) slot_alive[[node]] <<- FALSE
    }
    assign(key, state, envir = states)
    invisible(NULL)
  }

  collect_ready_results <- function(wait = FALSE) {
    collected <- FALSE
    busy_nodes <- which(slot_busy)
    if (length(busy_nodes) == 0L) return(FALSE)

    if (isTRUE(wait)) {
      try(slots[[busy_nodes[[1L]]]]$poll_process(poll_interval_ms), silent = TRUE)
    }
    for (node in busy_nodes) {
      ready <- tryCatch(
        identical(slots[[node]]$poll_process(0), "ready"),
        error = function(e) TRUE
      )
      if (!isTRUE(ready)) next
      result <- read_global_hal_slot_event(slots[[node]])
      if (is.null(result)) next
      process_slot_result(node, result)
      collected <- TRUE
    }
    collected
  }

  repeat {
    finalize_ready_runs()
    admit_runs()
    dispatch_jobs()
    finalize_ready_runs()

    done <- next_run > length(run_ids) &&
      length(active_keys) == 0L &&
      !any(slot_busy)
    if (isTRUE(done)) break

    if (any(slot_busy)) {
      if (!collect_ready_results(wait = FALSE)) {
        collect_ready_results(wait = TRUE)
      }
      next
    }

    # No task is running, yet work remains: mark active contexts as failed
    # instead of spinning forever (for example, after all slots crash).
    if (length(active_keys) > 0L) {
      for (key in active_keys) {
        state <- get(key, envir = states, inherits = FALSE)
        state$error <- "global scheduler deadlock: no live eligible HAL slot"
        state$pending_labels <- character()
        assign(key, state, envir = states)
      }
      finalize_ready_runs()
      next
    }
  }

  ordered_keys <- as.character(run_ids)
  finalized <- finalized[ordered_keys]
  trace <- if (length(trace_rows) == 0L) {
    data.frame()
  } else {
    do.call(rbind, trace_rows)
  }

  list(
    finalized = finalized,
    trace = trace,
    max_observed_global_fits = max_observed_global,
    max_observed_active_runs = max_observed_active,
    max_observed_per_run_fits = max_observed_per_run
  )
}
