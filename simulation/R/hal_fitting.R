###############################################################################
# HAL nuisance model fitting helpers
###############################################################################


#' Create a deterministic seed offset for a fit label.
#'
#' @param label Character label for a nuisance fit.
#'
#' @return Integer seed offset derived from the label.
stable_label_seed <- function(label) {
  chars <- utf8ToInt(label)
  as.integer(sum(chars * seq_along(chars)) %% 100000L)
}


#' Create reproducible HAL cross-validation fold IDs.
#'
#' The helper temporarily sets the random seed and then restores the previous
#' random-number state so fold construction does not perturb the rest of the
#' simulation replicate.
#'
#' @param n Integer number of observations.
#' @param nfolds Integer number of CV folds.
#' @param seed Integer seed used to shuffle fold assignments.
#'
#' @return Integer vector of fold IDs with length `n`.
make_foldid <- function(n, nfolds, seed) {
  n <- suppressWarnings(as.integer(n))
  nfolds <- suppressWarnings(as.integer(nfolds))
  if (is.na(n) || n < 3L || is.na(nfolds) ||
      nfolds < 3L || nfolds > n) {
    stop(
      "nfolds must be an integer between 3 and the number of observations.",
      call. = FALSE
    )
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  sample(rep(seq_len(nfolds), length.out = n))
}


#' Fit one HAL model and record timing.
#'
#' This is a wrapper around `hal9001::fit_hal()`. It does not implement HAL
#' itself; it manages CV workers, optional fit caching, and timing diagnostics
#' for one nuisance model.
#'
#' @param label Character label for the nuisance fit.
#' @param X Data frame or matrix of predictors.
#' @param Y Response vector.
#' @param weights Optional observation weights.
#' @param family Character outcome family passed to `hal9001::fit_hal()`.
#' @param smoothness_orders HAL smoothness order setting.
#' @param max_degree HAL maximum interaction degree.
#' @param num_knots HAL knot setting.
#' @param foldid Optional vector of CV fold IDs.
#' @param cv_folds Integer number of CV folds.
#' @param cv_nlambda Integer number of lambda values.
#' @param cv_workers Integer number of workers for HAL CV.
#' @param cv_cluster Optional reusable CV PSOCK cluster. It is only supported
#'   when complete HAL fits are running serially in the calling R process.
#' @param cache_file Optional `.rds` cache file path.
#' @param cache_key Expected v2 content key for `cache_file`.
#' @param use_cache Logical, whether to read/write the cache file.
#'
#' @return A list with `fit`, the fitted HAL object, and `timing`, a one-row data
#'   frame of fit diagnostics.
fit_hal_timed <- function(label, X, Y, weights = NULL, family = "binomial",
                          smoothness_orders = 1, max_degree = 2,
                          num_knots = REFERENCE_DEFAULT_NUM_KNOTS, foldid = NULL,
                          cv_folds = REFERENCE_DEFAULT_CV_FOLDS, cv_nlambda = 50,
                          cv_workers = 4, cv_cluster = NULL,
                          cache_file = NULL, cache_key = NULL,
                          use_cache = TRUE) {
  set_single_threaded_math()

  cv_folds <- suppressWarnings(as.integer(cv_folds))
  cv_nlambda <- suppressWarnings(as.integer(cv_nlambda))
  if (is.na(cv_folds) || cv_folds < 3L || cv_folds > nrow(X)) {
    stop(
      "cv_folds must be between 3 and the number of observations.",
      call. = FALSE
    )
  }
  if (is.na(cv_nlambda) || cv_nlambda < 1L) {
    stop("cv_nlambda must be a positive integer.", call. = FALSE)
  }

  cache_lookup_elapsed_sec <- 0
  if (isTRUE(use_cache) && !is.null(cache_file) && !is.null(cache_key)) {
    cache_t0 <- proc.time()
    cached <- read_fit_cache(cache_file, cache_key)
    cache_lookup_elapsed_sec <- unname(
      (proc.time() - cache_t0)[["elapsed"]]
    )
  } else {
    cached <- NULL
  }
  if (!is.null(cached)) {
    if (!"original_fit_elapsed_sec" %in% names(cached$timing)) {
      cached$timing$original_fit_elapsed_sec <- cached$timing$elapsed_sec
    }
    cached$timing$cache_lookup_elapsed_sec <- cache_lookup_elapsed_sec
    cached$timing$effective_elapsed_sec <- cache_lookup_elapsed_sec
    cached$timing$cache_read_pid <- Sys.getpid()
    cached$timing$cache_read_time <- format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
    cached$timing$cache_hit <- TRUE
    cached$timing$cache_file <- cache_file
    cached$timing$cache_key <- cache_key
    return(cached)
  }

  cv_workers <- min(max(1L, as.integer(cv_workers)), cv_folds)
  use_cv_parallel <- cv_workers > 1L
  owns_cv_cluster <- FALSE

  if (isTRUE(use_cv_parallel)) {
    cl <- cv_cluster
    if (is.null(cl)) {
      cl <- parallel::makeCluster(cv_workers)
      owns_cv_cluster <- TRUE
    } else if (length(cl) != cv_workers) {
      stop(
        "Reusable CV cluster size must equal cv_workers.",
        call. = FALSE
      )
    }
    doParallel::registerDoParallel(cl)
    on.exit({
      foreach::registerDoSEQ()
      if (isTRUE(owns_cv_cluster)) {
        try(parallel::stopCluster(cl), silent = TRUE)
      }
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  start_time <- Sys.time()
  t0 <- proc.time()
  fit <- hal9001::fit_hal(
    X = X,
    Y = Y,
    weights = weights,
    family = family,
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots,
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-4,
      nfolds = cv_folds,
      nlambda = cv_nlambda,
      foldid = foldid,
      parallel = use_cv_parallel,
      prediction_bounds = "default"
    )
  )
  elapsed_sec <- unname((proc.time() - t0)[["elapsed"]])
  end_time <- Sys.time()

  hal_times <- fit$times
  timing <- data.frame(
    fit_label = label,
    n_obs = nrow(X),
    n_covariates = ncol(X),
    response_mean = mean(Y),
    weight_sum = if (is.null(weights)) NA_real_ else sum(weights),
    elapsed_sec = elapsed_sec,
    elapsed_min = elapsed_sec / 60,
    original_fit_elapsed_sec = elapsed_sec,
    cache_lookup_elapsed_sec = cache_lookup_elapsed_sec,
    effective_elapsed_sec = elapsed_sec,
    hal_total_sec = unname(hal_times["total", "elapsed"]),
    hal_enumerate_basis_sec = unname(hal_times["enumerate_basis", "elapsed"]),
    hal_design_matrix_sec = unname(hal_times["design_matrix", "elapsed"]),
    hal_lasso_sec = unname(hal_times["lasso", "elapsed"]),
    lambda_star = fit$lambda_star,
    cv_folds = cv_folds,
    cv_nlambda = cv_nlambda,
    cv_workers = cv_workers,
    pid = Sys.getpid(),
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    cache_hit = FALSE,
    cache_read_pid = NA_integer_,
    cache_read_time = NA_character_,
    cache_file = if (is.null(cache_file)) NA_character_ else cache_file,
    cache_key = if (is.null(cache_key)) NA_character_ else cache_key,
    stringsAsFactors = FALSE
  )

  out <- list(
    fit = fit,
    timing = timing,
    cache_meta = list(
      schema = FIT_CACHE_SCHEMA,
      key = cache_key
    )
  )
  if (isTRUE(use_cache) && !is.null(cache_file) && !is.null(cache_key)) {
    atomic_save_rds(
      out,
      cache_file,
      compress = FALSE,
      overwrite = FALSE,
      verify = FALSE
    )
  }
  out
}


#' Run one named HAL nuisance-fit job.
#'
#' @param spec List describing one fit job, as produced by
#'   `make_reference_fit_jobs()`.
#' @param use_cache Logical, whether cached HAL fits may be used.
#' @param cv_cluster Optional reusable CV PSOCK cluster.
#'
#' @return Output from `fit_hal_timed()`.
run_one_hal_job <- function(spec, use_cache = TRUE, cv_cluster = NULL) {
  fit_hal_timed(
    label = spec$label,
    X = spec$X,
    Y = spec$Y,
    weights = spec$weights,
    family = spec$family,
    smoothness_orders = spec$smoothness_orders,
    max_degree = spec$max_degree,
    num_knots = spec$num_knots,
    foldid = spec$foldid,
    cv_folds = spec$cv_folds,
    cv_nlambda = spec$cv_nlambda,
    cv_workers = spec$cv_workers,
    cv_cluster = cv_cluster,
    cache_file = spec$cache_file,
    cache_key = spec$cache_key,
    use_cache = use_cache
  )
}


#' Resolve the complete-fit dispatch backend for this platform.
#'
#' @param backend Character backend: `"auto"`, `"fork"`, or `"psock"`.
#'
#' @return Resolved backend name.
resolve_fit_backend <- function(backend = "auto") {
  if (!backend %in% c("auto", "fork", "psock")) {
    stop("backend must be one of: auto, fork, psock", call. = FALSE)
  }
  if (identical(backend, "auto")) {
    return(if (.Platform$OS.type == "windows") "psock" else "fork")
  }
  if (identical(backend, "fork") && .Platform$OS.type == "windows") {
    stop("The fork backend is not available on Windows.", call. = FALSE)
  }
  backend
}


#' Create and initialize a reusable PSOCK cluster for complete HAL jobs.
#'
#' @param worker_count Positive integer number of workers.
#'
#' @return A PSOCK cluster initialized with packages, thread limits, and helper
#'   functions required by `run_one_hal_job()`.
make_hal_fit_psock_cluster <- function(worker_count) {
  worker_count <- max(1L, as.integer(worker_count))
  function_env <- environment(run_one_hal_job)
  package_name <- if (isNamespace(function_env)) {
    getNamespaceName(function_env)
  } else {
    NULL
  }
  package_mode <- identical(package_name, "popart")
  cl <- parallel::makeCluster(worker_count)
  initialized <- FALSE
  on.exit({
    if (!isTRUE(initialized)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)

  parallel::clusterCall(
    cl,
    function(package_mode, package_name) {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
      if (isTRUE(package_mode)) {
        loadNamespace(package_name)
      }
      NULL
    },
    package_mode = package_mode,
    package_name = package_name
  )
  if (isTRUE(package_mode)) {
    # Functions serialized from an installed package retain their namespace
    # environment. Loading that namespace in every worker avoids copying
    # implementation helpers into `.GlobalEnv`.
    parallel::clusterCall(
      cl,
      function(package_name) {
        loadNamespace("hal9001")
        loadNamespace("doParallel")
        loadNamespace(package_name)
        NULL
      },
      package_name = package_name
    )
  } else {
    # Preserve the standalone source-tree execution path.
    parallel::clusterEvalQ(cl, {
      loadNamespace("hal9001")
      loadNamespace("doParallel")
      NULL
    })

    helper_names <- c(
      "FIT_CACHE_SCHEMA",
      "fit_hal_timed",
      "run_one_hal_job",
      "set_single_threaded_math",
      "read_fit_cache",
      "quarantine_invalid_file",
      "safe_path_component",
      "recover_atomic_target",
      "atomic_promote_file",
      "atomic_save_rds"
    )
    parallel::clusterExport(
      cl,
      varlist = helper_names,
      envir = function_env
    )
  }

  attr(cl, "aipw_package_name") <- if (isTRUE(package_mode)) {
    package_name
  } else {
    NULL
  }

  initialized <- TRUE
  cl
}


#' Create a reusable PSOCK cluster for CV folds.
#'
#' This cluster is used only when complete HAL fits run serially in the main R
#' process. Socket cluster objects must not be sent to other PSOCK workers.
#' @noRd
make_hal_cv_psock_cluster <- function(worker_count) {
  worker_count <- max(1L, as.integer(worker_count))
  cl <- parallel::makeCluster(worker_count)
  initialized <- FALSE
  on.exit({
    if (!isTRUE(initialized)) {
      try(parallel::stopCluster(cl), silent = TRUE)
    }
  }, add = TRUE)
  parallel::clusterCall(
    cl,
    function() {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
      NULL
    }
  )
  initialized <- TRUE
  cl
}


#' Dispatch multiple HAL nuisance-fit jobs locally.
#'
#' Runs the nuisance fits serially, with forked workers, or with a PSOCK cluster
#' depending on `fit_workers` and `backend`.
#'
#' @param job_specs Named list of HAL fit specifications.
#' @param fit_workers Integer number of complete HAL fits to run concurrently.
#' @param backend Character dispatch backend: `"auto"`, `"fork"`, or `"psock"`.
#' @param use_cache Logical, whether cached HAL fits may be used.
#' @param fit_cluster Optional reusable outer PSOCK cluster.
#' @param cv_cluster Optional reusable CV PSOCK cluster. Only used when
#'   `fit_workers == 1`.
#'
#' @return Named list of outputs from `fit_hal_timed()`.
run_hal_jobs <- function(job_specs, fit_workers, backend = "auto",
                         use_cache = TRUE, fit_cluster = NULL,
                         cv_cluster = NULL) {
  n_jobs <- length(job_specs)
  worker_count <- min(max(1L, fit_workers), n_jobs)

  if (worker_count == 1L) {
    out <- lapply(
      job_specs,
      run_one_hal_job,
      use_cache = use_cache,
      cv_cluster = cv_cluster
    )
    names(out) <- names(job_specs)
    return(out)
  }

  backend <- resolve_fit_backend(backend)

  if (identical(backend, "fork")) {
    out <- parallel::mclapply(
      job_specs,
      run_one_hal_job,
      use_cache = use_cache,
      mc.cores = worker_count,
      mc.set.seed = FALSE,
      mc.preschedule = FALSE
    )
  } else if (identical(backend, "psock")) {
    owns_fit_cluster <- is.null(fit_cluster)
    cl <- if (isTRUE(owns_fit_cluster)) {
      make_hal_fit_psock_cluster(worker_count)
    } else {
      fit_cluster
    }
    if (length(cl) != worker_count) {
      stop(
        "Reusable fit cluster size must equal fit_workers.",
        call. = FALSE
      )
    }
    if (isTRUE(owns_fit_cluster)) {
      on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    }
    function_env <- environment(run_one_hal_job)
    package_name <- if (isNamespace(function_env)) {
      getNamespaceName(function_env)
    } else {
      NULL
    }
    package_mode <- identical(package_name, "popart")
    if (isTRUE(package_mode)) {
      # A caller may supply its own reusable PSOCK cluster, so make namespace
      # initialization explicit here as well as in the cluster constructor.
      if (!identical(attr(cl, "aipw_package_name"), package_name)) {
        parallel::clusterCall(
          cl,
          function(package_name) {
            loadNamespace(package_name)
            NULL
          },
          package_name = package_name
        )
      }
      out <- parallel::parLapplyLB(
        cl,
        job_specs,
        function(spec, use_cache, package_name) {
          run_job <- get(
            "run_one_hal_job",
            envir = asNamespace(package_name),
            inherits = FALSE
          )
          run_job(spec, use_cache = use_cache)
        },
        use_cache = use_cache,
        package_name = package_name
      )
    } else {
      out <- parallel::parLapplyLB(
        cl,
        job_specs,
        run_one_hal_job,
        use_cache = use_cache
      )
    }
  } else {
    stop("backend must be one of: auto, fork, psock", call. = FALSE)
  }

  names(out) <- names(job_specs)
  out
}
