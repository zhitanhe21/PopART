# Global HAL scheduler --------------------------------------------------------

# One fixed worker pool fits all HAL jobs from the active Monte Carlo
# replicates. Replicates are prepared in batches to limit memory use.

run_global_scheduler <- function(
    run_ids,
    prepare_run,
    finalize_run,
    global_fit_slots,
    active_replicates) {

  # Start one worker pool ---------------------------------------------------

  fit_job <- .fit_hal_job
  environment(fit_job) <- baseenv()
  cluster <- parallel::makePSOCKcluster(global_fit_slots)
  on.exit(parallel::stopCluster(cluster))

  parallel::clusterCall(cluster, function() {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1"
    )
    NULL
  })

  # Prepare, fit, and finish each batch ------------------------------------

  batches <- split(
    run_ids,
    ceiling(seq_along(run_ids) / active_replicates)
  )
  results <- list()

  for (batch in batches) {
    contexts <- lapply(batch, prepare_run)
    job_index <- do.call(rbind, lapply(seq_along(contexts), function(i) {
      data.frame(
        replicate = i,
        fit = names(contexts[[i]]$jobs),
        stringsAsFactors = FALSE
      )
    }))
    jobs <- unlist(lapply(contexts, `[[`, "jobs"), recursive = FALSE)
    fits <- parallel::parLapplyLB(cluster, jobs, fit_job)

    for (i in seq_along(contexts)) {
      replicate_fits <- fits[job_index$replicate == i]
      names(replicate_fits) <- job_index$fit[job_index$replicate == i]
      context <- contexts[[i]]
      context$jobs <- NULL
      run_id <- batch[[i]]
      results[[as.character(run_id)]] <- finalize_run(
        run_id, context, replicate_fits
      )
    }
  }

  list(results = results[as.character(run_ids)])
}
