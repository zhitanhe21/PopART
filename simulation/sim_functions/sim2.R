# Simulation 2: global-scheduler HAL AIPW -----------------------------------

# 1. Generate one trial and one auxiliary sample.
generate_sim2_data <- function(
    n_clusters, n_trial, n_auxiliary, p_response, p_censoring, seed) {

  n_trial <- n_clusters * round(n_trial / n_clusters)
  n_auxiliary <- n_clusters * round(n_auxiliary / n_clusters)
  set.seed(seed)

  cluster_risk <- c(
    0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
    1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76
  )
  treatment_by_cluster <- sample(rep(c(0, 1), n_clusters / 2))
  response_intercept <- 0.4449358 + qlogis(p_response) - qlogis(0.5)
  censoring_intercept <- -0.8711302 + qlogis(p_censoring) - qlogis(0.3)

  cluster <- rep(seq_len(n_clusters), length.out = n_trial)
  A <- treatment_by_cluster[cluster]
  X1 <- cluster_risk[cluster]
  W1 <- rbinom(n_trial, 1, 0.5)
  W2 <- rnorm(n_trial)
  W3 <- rnorm(n_trial)
  R <- rbinom(
    n_trial, 1,
    plogis(response_intercept - 2 * A * W1 + 0.2 * W2 * W3^2)
  )
  C <- rbinom(
    n_trial, 1,
    plogis(censoring_intercept - 0.25 * A + 0.25 * W1 + 0.2 * W2^2 * W3)
  )
  C[R == 0L] <- 0L
  mu0 <- plogis(-1 + 2 * W1)
  mu1 <- plogis(-W1 + 0.2 * sin(W2) - 0.2 * W3^2 + 0.25 * X1)
  Y0 <- rbinom(n_trial, 1, mu0)
  Y1 <- rbinom(n_trial, 1, mu1)
  Y <- (1 - A) * Y0 + A * Y1
  Y[R == 0L | C == 1L] <- 0L

  trial <- data.frame(
    id = seq_len(n_trial), cluster, X1, W1, W2, W3, A, R, C, Y, wt = 1,
    S = 1L
  )

  cluster <- rep(seq_len(n_clusters), length.out = n_auxiliary)
  W1 <- rbinom(n_auxiliary, 1, 0.75)
  auxiliary <- data.frame(
    id = n_trial + seq_len(n_auxiliary),
    cluster,
    X1 = cluster_risk[cluster],
    W1,
    W2 = rnorm(n_auxiliary),
    W3 = rnorm(n_auxiliary),
    A = 0L,
    R = 0L,
    C = 0L,
    Y = 0L,
    wt = ifelse(W1 == 0L, 0.75, 0.25),
    S = 0L
  )
  auxiliary$wt <- auxiliary$wt / mean(auxiliary$wt)

  list(
    data = rbind(trial, auxiliary),
    truth = c(eta0 = mean(mu0), eta1 = mean(mu1)),
    n_trial = n_trial,
    n_auxiliary = n_auxiliary
  )
}


# 2. Run the HAL AIPW study through the global scheduler.
run_sim2 <- function(
    sample_sizes,
    run_id,
    out_dir,
    n_clusters = 20L,
    p_response = 0.5,
    p_censoring = 0.3,
    mc_reps = 20L,
    base_seed = 91000000L,
    global_fit_slots = 4L,
    active_replicates = 2L,
    cv_folds = 5L,
    cv_lambdas = 50L) {

  grid <- sample_sizes[rep(seq_len(nrow(sample_sizes)), mc_reps), ]
  grid$replicate <- rep(seq_len(mc_reps), each = nrow(sample_sizes))
  grid$seed <- base_seed + seq_len(nrow(grid))
  rownames(grid) <- NULL

  # Scheduler callbacks ----------------------------------------------------

  prepare <- function(i) {
    generated <- generate_sim2_data(
      n_clusters, grid$n_trial[[i]], grid$n_auxiliary[[i]],
      p_response, p_censoring, grid$seed[[i]]
    )
    covariates <- c("X1", "W1", "W2", "W3")
    list(
      jobs = .make_popart_jobs(
        generated$data, covariates, cv_folds, cv_lambdas, grid$seed[[i]],
        generated$data$cluster
      ),
      data = generated$data,
      covariates = covariates,
      truth = generated$truth,
      seed = grid$seed[[i]],
      n_trial = generated$n_trial,
      n_auxiliary = generated$n_auxiliary
    )
  }

  finish <- function(i, context, fits) {
    popart <- .finish_popart(
      context$data, context$covariates, fits,
      cluster = context$data$cluster
    )
    do.call(rbind, lapply(c("Naive", "Proposed"), function(version) {
      estimate <- popart$estimates[
        popart$estimates$estimator == "AIPW" &
          popart$estimates$version == version,
      ]
      variance <- popart$variance[
        popart$variance$estimator == "AIPW" &
          popart$variance$version == version,
      ]
      covariance <- popart$covariance[[paste0("aipw_", tolower(version))]]
      eta0 <- estimate$estimate[estimate$parameter == "mean_control"]
      eta1 <- estimate$estimate[estimate$parameter == "mean_treated"]
      data.frame(
        Estimator = "AIPW", Version = version,
        etahat_0 = eta0, etahat_1 = eta1,
        rdhat = eta1 - eta0, rrhat = eta1 / eta0,
        cov_00 = variance$variance[variance$parameter == "mean_control"],
        cov_01 = covariance[1, 2],
        cov_11 = variance$variance[variance$parameter == "mean_treated"],
        var_rd = variance$variance[variance$parameter == "risk_difference"],
        var_rr = variance$variance[variance$parameter == "risk_ratio"],
        eta_0 = context$truth[["eta0"]],
        eta_1 = context$truth[["eta1"]],
        rd = context$truth[["eta1"]] - context$truth[["eta0"]],
        rr = context$truth[["eta1"]] / context$truth[["eta0"]],
        seed = context$seed, n_trial = context$n_trial,
        n_aux = context$n_auxiliary,
        run_id = run_id, replicate = grid$replicate[[i]]
      )
    }))
  }

  # Global HAL fitting -----------------------------------------------------

  cat(
    "Running", mc_reps, "Monte Carlo replicates for each of",
    nrow(sample_sizes), "sample-size combinations\n"
  )
  scheduler <- run_global_scheduler(
    run_ids = seq_len(nrow(grid)),
    prepare_run = prepare,
    finalize_run = finish,
    global_fit_slots = global_fit_slots,
    active_replicates = active_replicates
  )
  results <- do.call(rbind, unname(scheduler$results))
  rownames(results) <- NULL

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(out_dir, paste0("monte_carlo_", run_id, "_results.csv"))
  write.csv(results, result_file, row.names = FALSE)
  cat("Results:", result_file, "\n")

  invisible(list(results = results, result_file = result_file))
}
