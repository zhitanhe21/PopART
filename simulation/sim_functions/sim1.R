# Simulation 1: parametric model misspecification ---------------------------

# 1. Generate one trial and one auxiliary sample.
generate_sim1_data <- function(
    n_clusters, n_trial, n_auxiliary, p_response, p_censoring, seed) {

  response_intercept <- stats::uniroot(
    function(value) {
      mean(stats::plogis(c(value, value - 2, value, value))) - p_response
    },
    c(-20, 20)
  )$root
  combinations <- expand.grid(A = c(0, 1), W1 = c(0, 1))
  response_probability <- stats::plogis(
    response_intercept + combinations$A - combinations$W1
  )
  censoring_intercept <- stats::uniroot(
    function(value) {
      censoring_probability <- stats::plogis(
        value - 0.25 * combinations$A + 0.25 * combinations$W1
      )
      sum(censoring_probability * response_probability * 0.25) /
        sum(response_probability * 0.25) - p_censoring
    },
    c(-20, 20)
  )$root

  n_trial <- n_clusters * round(n_trial / n_clusters)
  n_auxiliary <- n_clusters * round(n_auxiliary / n_clusters)
  set.seed(seed)
  cluster_risk <- c(
    0.92, 0.78, 0.07, -1.99, 0.62, -0.06, -0.16, -1.47, -0.48, 0.42,
    1.36, -0.10, 0.39, -0.05, -1.38, -0.41, -0.39, -0.06, 1.10, 0.76
  )
  treatment_by_cluster <- sample(rep(c(0L, 1L), n_clusters / 2L))

  cluster <- rep(seq_len(n_clusters), length.out = n_trial)
  A <- treatment_by_cluster[cluster]
  X1 <- cluster_risk[cluster]
  W1 <- stats::rbinom(n_trial, 1, 0.5)
  W2 <- stats::rnorm(n_trial)
  R <- stats::rbinom(
    n_trial, 1, stats::plogis(response_intercept - 2 * A * W1)
  )
  C <- stats::rbinom(
    n_trial, 1, stats::plogis(censoring_intercept - 0.25 * A + 0.25 * W1)
  )
  C[R == 0L] <- 0L
  probability0 <- stats::plogis(-1 + 2 * W1 + 0.5 * W2 + 0.25 * X1)
  probability1 <- stats::plogis(-W1 - 0.5 * W2)
  Y0 <- stats::rbinom(n_trial, 1, probability0)
  Y1 <- stats::rbinom(n_trial, 1, probability1)
  Y <- (1 - A) * Y0 + A * Y1
  Y[R == 0L | C == 1L] <- 0L
  trial <- data.frame(
    id = seq_len(n_trial), cluster, X1, W1, W2, A, R, C, Y, wt = 1, S = 1L
  )

  cluster <- rep(seq_len(n_clusters), length.out = n_auxiliary)
  A <- treatment_by_cluster[cluster]
  W1 <- stats::rbinom(n_auxiliary, 1, 0.75)
  auxiliary <- data.frame(
    id = n_trial + seq_len(n_auxiliary),
    cluster,
    X1 = cluster_risk[cluster],
    W1,
    W2 = stats::rnorm(n_auxiliary),
    A,
    R = 0L,
    C = 0L,
    Y = 0L,
    wt = ifelse(W1 == 0L, 0.75, 0.25),
    S = 0L
  )

  list(
    data = rbind(trial, auxiliary),
    truth = c(eta0 = mean(Y0), eta1 = mean(Y1)),
    n_trial = n_trial,
    n_auxiliary = n_auxiliary
  )
}
# 2. Run the Monte Carlo study and fit the six parametric estimators.
run_sim1 <- function(
    sample_sizes, run_id, out_dir,
    n_clusters = 20L, p_response = 0.5, p_censoring = 0.3,
    mc_reps = 30L, base_seed = 11000000L) {

  # Monte Carlo settings ---------------------------------------------------

  grid <- sample_sizes[rep(seq_len(nrow(sample_sizes)), mc_reps), ]
  grid$replicate <- rep(seq_len(mc_reps), each = nrow(sample_sizes))
  grid$seed <- base_seed + seq_len(nrow(grid))
  rownames(grid) <- NULL
  scenarios <- expand.grid(
    mu_correct = c(TRUE, FALSE),
    pi_correct = c(TRUE, FALSE)
  )

  cat(
    "Running", mc_reps, "Simulation 1 replicates for each of",
    nrow(sample_sizes), "sample-size combinations\n"
  )

  # Simulation pipeline ----------------------------------------------------

  results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    generated <- generate_sim1_data(
      n_clusters, grid$n_trial[[i]], grid$n_auxiliary[[i]],
      p_response, p_censoring, grid$seed[[i]]
    )

    do.call(rbind, lapply(seq_len(nrow(scenarios)), function(j) {
      mu_correct <- scenarios$mu_correct[[j]]
      pi_correct <- scenarios$pi_correct[[j]]
      outcome_formula <- if (mu_correct) {
        Y ~ A * (X1 + W1 + W2)
      } else {
        Y ~ A * (X1 + W2)
      }
      propensity_formula <- if (pi_correct) {
        Q ~ X1 + W1 + W2
      } else {
        Q ~ X1 + W2
      }
      censoring_formula <- if (pi_correct) {
        C ~ X1 + W1 + W2
      } else {
        C ~ X1 + W2
      }

      estimates <- rbind(
        fit_gformula(generated$data, outcome_formula, "naive"),
        fit_gformula(generated$data, outcome_formula, "proposed"),
        fit_ipw(generated$data, censoring_formula, "naive"),
        fit_ipw(generated$data, propensity_formula, "proposed"),
        fit_aipw(generated$data, outcome_formula, censoring_formula, "naive"),
        fit_aipw(generated$data, outcome_formula, propensity_formula, "proposed")
      )
      estimates$Estimator <- rep(c("G-Formula", "IPW", "AIPW"), each = 2L)
      estimates$Version <- rep(c("Naive", "Proposed"), 3L)
      estimates$rdhat <- estimates$etahat_1 - estimates$etahat_0
      estimates$rrhat <- estimates$etahat_1 / estimates$etahat_0
      estimates$var_rd <- estimates$cov_11 + estimates$cov_00 -
        2 * estimates$cov_01
      estimates$var_rr <- estimates$cov_11 / estimates$etahat_0^2 +
        estimates$cov_00 * estimates$etahat_1^2 / estimates$etahat_0^4 -
        2 * estimates$cov_01 * estimates$etahat_1 / estimates$etahat_0^3
      estimates <- estimates[c(
        "Estimator", "Version", "etahat_0", "etahat_1", "rdhat", "rrhat",
        "cov_00", "cov_01", "cov_11", "var_rd", "var_rr"
      )]

      transform(
        estimates,
        eta_0 = generated$truth[["eta0"]],
        eta_1 = generated$truth[["eta1"]],
        rd = generated$truth[["eta1"]] - generated$truth[["eta0"]],
        rr = generated$truth[["eta1"]] / generated$truth[["eta0"]],
        seed = grid$seed[[i]],
        n_trial = generated$n_trial,
        n_aux = generated$n_auxiliary,
        mu_correct = mu_correct,
        pi_correct = pi_correct,
        run_id = run_id,
        replicate = grid$replicate[[i]]
      )
    }))
  }))

  # Save results -----------------------------------------------------------

  rownames(results) <- NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(out_dir, paste0("simulation1_", run_id, "_results.csv"))
  utils::write.csv(results, result_file, row.names = FALSE)
  invisible(list(results = results, result_file = result_file))
}
