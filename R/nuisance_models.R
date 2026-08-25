###############################################################################
# Nuisance-model fitting
###############################################################################


.make_cv_fold_ids <- function(n, n_folds, seed) {
  n <- suppressWarnings(as.integer(n))
  n_folds <- suppressWarnings(as.integer(n_folds))
  seed <- suppressWarnings(as.integer(seed))
  if (length(n) != 1L || is.na(n) || n < 3L ||
      length(n_folds) != 1L || is.na(n_folds) ||
      n_folds < 3L || n_folds > n) {
    stop(
      "n_folds must be an integer between 3 and the number of observations.",
      call. = FALSE
    )
  }
  if (length(seed) != 1L || is.na(seed) || seed < 0L) {
    stop("seed must be a nonnegative integer.", call. = FALSE)
  }

  previous_seed <- if (exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(previous_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  sample(rep(seq_len(n_folds), length.out = n))
}


.limit_math_threads <- function() {
  variables <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS"
  )
  previous <- Sys.getenv(variables, unset = NA_character_)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  previous
}


.restore_math_threads <- function(previous) {
  missing <- is.na(previous)
  if (any(missing)) {
    Sys.unsetenv(names(previous)[missing])
  }
  if (any(!missing)) {
    do.call(Sys.setenv, as.list(previous[!missing]))
  }
  invisible(TRUE)
}


.make_cv_cluster <- function(n_workers) {
  n_workers <- max(1L, as.integer(n_workers))
  cluster <- parallel::makeCluster(n_workers)
  initialized <- FALSE
  on.exit({
    if (!initialized) {
      try(parallel::stopCluster(cluster), silent = TRUE)
    }
  }, add = TRUE)
  parallel::clusterCall(cluster, function() {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    NULL
  })
  initialized <- TRUE
  cluster
}


.hal_time_value <- function(times, row) {
  if (is.null(times) || !is.matrix(times) ||
      !row %in% rownames(times) || !"elapsed" %in% colnames(times)) {
    return(NA_real_)
  }
  unname(times[row, "elapsed"])
}


.fit_hal_job <- function(spec, cv_cluster = NULL) {
  previous_threads <- .limit_math_threads()
  on.exit(.restore_math_threads(previous_threads), add = TRUE)

  n_cv_folds <- as.integer(spec$n_cv_folds)
  n_lambda_values <- as.integer(spec$n_lambda_values)
  n_cv_workers <- min(
    max(1L, as.integer(spec$n_cv_workers)),
    n_cv_folds
  )
  use_parallel_cv <- n_cv_workers > 1L
  owns_cluster <- FALSE

  if (use_parallel_cv) {
    cluster <- cv_cluster
    if (is.null(cluster)) {
      cluster <- .make_cv_cluster(n_cv_workers)
      owns_cluster <- TRUE
    } else if (length(cluster) != n_cv_workers) {
      stop(
        "The reusable CV cluster size must equal n_cv_workers.",
        call. = FALSE
      )
    }
    doParallel::registerDoParallel(cluster)
    on.exit({
      foreach::registerDoSEQ()
      if (owns_cluster) {
        try(parallel::stopCluster(cluster), silent = TRUE)
      }
    }, add = TRUE)
  } else {
    foreach::registerDoSEQ()
    on.exit(foreach::registerDoSEQ(), add = TRUE)
  }

  started_at <- Sys.time()
  timer <- proc.time()
  fit <- hal9001::fit_hal(
    X = spec$X,
    Y = spec$Y,
    weights = spec$weights,
    family = "binomial",
    smoothness_orders = spec$smoothness_orders,
    max_degree = spec$max_degree,
    num_knots = spec$num_knots,
    return_lasso = isTRUE(spec$return_lasso),
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-4,
      nfolds = n_cv_folds,
      nlambda = n_lambda_values,
      foldid = spec$fold_ids,
      parallel = use_parallel_cv,
      prediction_bounds = "default"
    )
  )
  elapsed_seconds <- unname((proc.time() - timer)[["elapsed"]])
  finished_at <- Sys.time()

  list(
    fit = fit,
    timing = data.frame(
      outer_fold = as.integer(spec$outer_fold),
      fit = spec$label,
      learner = "hal",
      observations = nrow(spec$X),
      predictors = ncol(spec$X),
      response_mean = mean(spec$Y),
      weight_sum = if (is.null(spec$weights)) NA_real_ else sum(spec$weights),
      elapsed_seconds = elapsed_seconds,
      basis_seconds = .hal_time_value(fit$times, "enumerate_basis"),
      design_matrix_seconds = .hal_time_value(fit$times, "design_matrix"),
      lasso_seconds = .hal_time_value(fit$times, "lasso"),
      selected_lambda = unname(fit$lambda_star),
      n_cv_folds = n_cv_folds,
      n_lambda_values = n_lambda_values,
      n_cv_workers = n_cv_workers,
      n_boosting_rounds = NA_integer_,
      xgboost_max_depth = NA_integer_,
      xgboost_learning_rate = NA_real_,
      process_id = Sys.getpid(),
      started_at = format(started_at, "%Y-%m-%d %H:%M:%S"),
      finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  )
}


.fit_xgboost_job <- function(spec) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop(
      "Package 'xgboost' is required when a nuisance model uses xgboost.",
      call. = FALSE
    )
  }
  previous_threads <- .limit_math_threads()
  on.exit(.restore_math_threads(previous_threads), add = TRUE)

  matrix_x <- as.matrix(spec$X)
  storage.mode(matrix_x) <- "double"
  dmatrix_arguments <- list(data = matrix_x, label = spec$Y)
  if (!is.null(spec$weights)) {
    dmatrix_arguments$weight <- spec$weights
  }
  training_data <- do.call(xgboost::xgb.DMatrix, dmatrix_arguments)

  started_at <- Sys.time()
  timer <- proc.time()
  fit <- xgboost::xgb.train(
    params = list(
      objective = "binary:logistic",
      booster = "gbtree",
      max_depth = as.integer(spec$xgboost_max_depth),
      eta = as.numeric(spec$xgboost_learning_rate),
      min_child_weight = 1,
      subsample = 1,
      colsample_bytree = 1,
      lambda = 1,
      tree_method = "hist",
      nthread = 1L,
      seed = as.integer(spec$xgboost_seed),
      verbosity = 0L
    ),
    data = training_data,
    nrounds = as.integer(spec$xgboost_nrounds),
    verbose = 0L
  )
  elapsed_seconds <- unname((proc.time() - timer)[["elapsed"]])
  finished_at <- Sys.time()

  list(
    fit = fit,
    timing = data.frame(
      outer_fold = as.integer(spec$outer_fold),
      fit = spec$label,
      learner = "xgboost",
      observations = nrow(spec$X),
      predictors = ncol(spec$X),
      response_mean = mean(spec$Y),
      weight_sum = if (is.null(spec$weights)) NA_real_ else sum(spec$weights),
      elapsed_seconds = elapsed_seconds,
      basis_seconds = NA_real_,
      design_matrix_seconds = NA_real_,
      lasso_seconds = NA_real_,
      selected_lambda = NA_real_,
      n_cv_folds = NA_integer_,
      n_lambda_values = NA_integer_,
      n_cv_workers = NA_integer_,
      n_boosting_rounds = as.integer(spec$xgboost_nrounds),
      xgboost_max_depth = as.integer(spec$xgboost_max_depth),
      xgboost_learning_rate = as.numeric(spec$xgboost_learning_rate),
      process_id = Sys.getpid(),
      started_at = format(started_at, "%Y-%m-%d %H:%M:%S"),
      finished_at = format(finished_at, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  )
}


.fit_nuisance_job <- function(spec, cv_cluster = NULL) {
  if (identical(spec$learner, "hal")) {
    return(.fit_hal_job(spec, cv_cluster = cv_cluster))
  }
  if (identical(spec$learner, "xgboost")) {
    return(.fit_xgboost_job(spec))
  }
  stop("Unsupported nuisance learner: ", spec$learner, ".", call. = FALSE)
}


.predict_nuisance_fit <- function(fit, new_data) {
  if (inherits(fit, "hal9001")) {
    return(as.numeric(stats::predict(fit, new_data = new_data)))
  }
  if (inherits(fit, "xgb.Booster")) {
    matrix_x <- as.matrix(new_data)
    storage.mode(matrix_x) <- "double"
    return(as.numeric(stats::predict(fit, newdata = matrix_x)))
  }
  as.numeric(stats::predict(fit, new_data = new_data))
}


.resolve_parallel_backend <- function(backend) {
  if (!backend %in% c("auto", "fork", "psock")) {
    stop(
      "parallel_backend must be one of: auto, fork, psock.",
      call. = FALSE
    )
  }
  if (identical(backend, "auto")) {
    return(if (.Platform$OS.type == "windows") "psock" else "fork")
  }
  if (identical(backend, "fork") && .Platform$OS.type == "windows") {
    stop("The fork backend is not available on Windows.", call. = FALSE)
  }
  backend
}


.make_fit_cluster <- function(n_workers) {
  cluster <- parallel::makeCluster(n_workers)
  initialized <- FALSE
  on.exit({
    if (!initialized) {
      try(parallel::stopCluster(cluster), silent = TRUE)
    }
  }, add = TRUE)

  controller_library_paths <- .libPaths()
  parallel::clusterCall(cluster, function(controller_library_paths) {
    .libPaths(unique(c(controller_library_paths, .libPaths())))
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    loadNamespace("popart")
    NULL
  }, controller_library_paths = controller_library_paths)
  initialized <- TRUE
  cluster
}


.run_nuisance_jobs <- function(job_specs, n_fit_workers,
                               parallel_backend = "auto",
                               cv_cluster = NULL) {
  if (!is.list(job_specs) || length(job_specs) < 1L || is.null(names(job_specs))) {
    stop("job_specs must be a nonempty named list.", call. = FALSE)
  }
  n_workers <- min(max(1L, as.integer(n_fit_workers)), length(job_specs))

  if (n_workers == 1L) {
    output <- lapply(job_specs, .fit_nuisance_job, cv_cluster = cv_cluster)
    names(output) <- names(job_specs)
    return(output)
  }
  if (!is.null(cv_cluster)) {
    stop(
      "A CV cluster cannot be supplied while complete fits run concurrently.",
      call. = FALSE
    )
  }

  backend <- .resolve_parallel_backend(parallel_backend)
  if (identical(backend, "fork")) {
    output <- parallel::mclapply(
      job_specs,
      .fit_nuisance_job,
      mc.cores = n_workers,
      mc.set.seed = FALSE,
      mc.preschedule = FALSE
    )
  } else {
    cluster <- .make_fit_cluster(n_workers)
    on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
    output <- parallel::parLapplyLB(
      cluster,
      job_specs,
      function(spec) {
        fitter <- get(
          ".fit_nuisance_job",
          envir = asNamespace("popart"),
          inherits = FALSE
        )
        fitter(spec)
      }
    )
  }

  names(output) <- names(job_specs)
  output
}
