###############################################################################
# AIPW nuisance tasks and estimator calculations
###############################################################################


#' Construct nuisance-model jobs for the reference estimator.
#'
#' Builds the four HAL tasks needed by the reference study: two proposed Q
#' regressions, one shared outcome regression reused by both AIPW versions, and
#' one censoring regression for the naive estimator.
#'
#' @param dat Combined reference-study data frame.
#' @param seed Integer simulation seed.
#' @param n_trial Integer rounded trial sample size.
#' @param n_aux Integer rounded auxiliary sample size.
#' @param config List returned by `configure_fitcv_parallel()`.
#' @param cache_dir Optional cache directory for HAL fit objects.
#'
#' @return A list containing analysis data, covariate names, and a named `jobs`
#'   list of HAL fit specifications.
make_reference_fit_jobs <- function(dat, seed, n_trial, n_aux, config,
                                    cache_dir = NULL) {
  mu_covariates <- c("X1", "W1", "W2", "W3")
  pi_covariates <- c("X1", "W1", "W2", "W3")
  mu_cols <- c("A", mu_covariates)

  dat_naive <- dat[dat$S == 1 & dat$R == 1, , drop = FALSE]
  dat$Q <- dat$S * dat$R * (1 - dat$C)

  idx_shared_outcome <- dat$S == 1 & dat$R == 1 & dat$C == 0
  idx_q0 <- (dat$Q == 1 & dat$A == 0) | dat$S == 0
  idx_q1 <- (dat$Q == 1 & dat$A == 1) | dat$S == 0

  make_fit_job <- function(label, X, Y, weights = NULL) {
    family <- "binomial"
    smoothness_orders <- 1
    max_degree <- 2
    # Experimental laptop default. This deliberately uses a coarser HAL
    # dictionary than the manuscript reference setting and is tracked in cache and
    # checkpoint identities.
    num_knots <- REFERENCE_DEFAULT_NUM_KNOTS
    fold_seed <- seed + stable_label_seed(label)
    foldid <- make_foldid(length(Y), config$cv_folds, fold_seed)
    cache_identity <- NULL
    cache_file <- NULL
    if (!is.null(cache_dir) && nzchar(cache_dir)) {
      cache_identity <- make_fit_cache_identity(
        label = label,
        X = X,
        Y = Y,
        weights = weights,
        family = family,
        smoothness_orders = smoothness_orders,
        max_degree = max_degree,
        num_knots = num_knots,
        foldid = foldid,
        cv_folds = config$cv_folds,
        cv_nlambda = config$cv_nlambda
      )
      cache_file <- make_fit_cache_file(
        cache_dir = cache_dir,
        label = label,
        seed = seed,
        n_trial = n_trial,
        n_aux = n_aux,
        cache_key = cache_identity$key
      )
    }

    list(
      label = label,
      X = X,
      Y = Y,
      weights = weights,
      family = family,
      smoothness_orders = smoothness_orders,
      max_degree = max_degree,
      num_knots = num_knots,
      foldid = foldid,
      cv_folds = config$cv_folds,
      cv_nlambda = config$cv_nlambda,
      cv_workers = config$cv_workers_per_fit,
      cache_file = cache_file,
      cache_key = if (is.null(cache_identity)) NULL else cache_identity$key
    )
  }

  list(
    dat_naive = dat_naive,
    dat_all = dat,
    mu_covariates = mu_covariates,
    pi_covariates = pi_covariates,
    jobs = list(
      proposed_Q_reg_0 = make_fit_job(
        "proposed_Q_reg_0",
        dat[idx_q0, pi_covariates, drop = FALSE],
        dat$Q[idx_q0],
        weights = dat$wt[idx_q0]
      ),
      proposed_Q_reg_1 = make_fit_job(
        "proposed_Q_reg_1",
        dat[idx_q1, pi_covariates, drop = FALSE],
        dat$Q[idx_q1],
        weights = dat$wt[idx_q1]
      ),
      shared_outcome_reg = make_fit_job(
        "shared_outcome_reg",
        dat[idx_shared_outcome, mu_cols, drop = FALSE],
        dat$Y[idx_shared_outcome]
      ),
      naive_censor_reg = make_fit_job(
        "naive_censor_reg",
        dat_naive[, pi_covariates, drop = FALSE],
        dat_naive$C
      )
    )
  )
}


#' Compute the naive AIPW estimates.
#'
#' Uses the shared outcome regression and the naive censoring regression to
#' estimate eta(0), eta(1), risk difference, risk ratio, and covariance terms.
#'
#' @param dat Trial responder data used by the naive estimator.
#' @param outcome_reg Fitted HAL outcome model.
#' @param censor_reg Fitted HAL censoring model.
#' @param mu_covariates Character vector of outcome-regression covariates.
#' @param pi_covariates Character vector of censoring-regression covariates.
#' @param pA Numeric probability of treatment assignment.
#'
#' @return One-row data frame with eta estimates and covariance estimates.
compute_naive_eta <- function(dat, outcome_reg, censor_reg, mu_covariates,
                              pi_covariates, pA = 0.5) {
  mu_cols <- c("A", mu_covariates)
  pred_idx <- dat$S == 0 | dat$R == 1

  new0 <- dat[pred_idx, mu_cols, drop = FALSE]
  new1 <- new0
  new0$A <- 0
  new1$A <- 1

  dat$muhat_0 <- NA_real_
  dat$muhat_1 <- NA_real_
  dat$muhat_0[pred_idx] <- predict(outcome_reg, new_data = new0)
  dat$muhat_1[pred_idx] <- predict(outcome_reg, new_data = new1)

  dat$piC <- 1 - predict(censor_reg, new_data = dat[, pi_covariates, drop = FALSE])
  dat$piA <- ifelse(dat$A == 1, pA, 1 - pA)
  dat$pihat <- dat$piA * dat$piC

  uncens_idx <- dat$C == 0
  n_trial_hat <- c(
    sum((1 - dat$A[uncens_idx]) / dat$pihat[uncens_idx]),
    sum(dat$A[uncens_idx] / dat$pihat[uncens_idx])
  )

  n_aux <- sum(dat$wt)
  n_dat <- nrow(dat)

  ipw_0 <- ifelse(
    dat$C == 0 & dat$A == 0,
    (dat$Y - dat$muhat_0) / dat$pihat,
    0
  )
  ipw_1 <- ifelse(
    dat$C == 0 & dat$A == 1,
    (dat$Y - dat$muhat_1) / dat$pihat,
    0
  )

  if_0 <- ((1 - dat$C) * ipw_0 / n_trial_hat[1] +
             dat$muhat_0 * dat$wt / n_aux) * n_dat
  if_1 <- ((1 - dat$C) * ipw_1 / n_trial_hat[2] +
             dat$muhat_1 * dat$wt / n_aux) * n_dat

  data.frame(
    etahat_0 = mean(if_0),
    etahat_1 = mean(if_1),
    cov_00 = var(if_0) / n_dat,
    cov_01 = cov(if_0, if_1) / n_dat,
    cov_11 = var(if_1) / n_dat
  )
}


#' Compute the proposed AIPW estimates.
#'
#' Uses the shared outcome regression and the two proposed Q regressions to
#' estimate eta(0), eta(1), risk difference, risk ratio, and covariance terms.
#'
#' @param dat Combined trial and auxiliary data.
#' @param outcome_reg Fitted HAL outcome model.
#' @param Q_reg_0 Fitted proposed Q model among arm 0 and auxiliary records.
#' @param Q_reg_1 Fitted proposed Q model among arm 1 and auxiliary records.
#' @param mu_covariates Character vector of outcome-regression covariates.
#' @param pi_covariates Character vector of Q-regression covariates.
#'
#' @return One-row data frame with eta estimates and covariance estimates.
compute_proposed_eta <- function(dat, outcome_reg, Q_reg_0, Q_reg_1,
                                 mu_covariates, pi_covariates) {
  mu_cols <- c("A", mu_covariates)
  dat$Q <- dat$S * dat$R * (1 - dat$C)

  pred_idx <- dat$S == 0 | dat$R == 1
  new0 <- dat[pred_idx, mu_cols, drop = FALSE]
  new1 <- new0
  new0$A <- 0
  new1$A <- 1

  dat$muhat_0 <- NA_real_
  dat$muhat_1 <- NA_real_
  dat$muhat_0[pred_idx] <- predict(outcome_reg, new_data = new0)
  dat$muhat_1[pred_idx] <- predict(outcome_reg, new_data = new1)

  Q_ind <- dat$Q == 1
  pred_q0_idx <- Q_ind & dat$A == 0
  pred_q1_idx <- Q_ind & dat$A == 1

  dat$Q_prob <- NA_real_
  if (any(pred_q0_idx)) {
    dat$Q_prob[pred_q0_idx] <- predict(
      Q_reg_0,
      new_data = dat[pred_q0_idx, pi_covariates, drop = FALSE]
    )
  }
  if (any(pred_q1_idx)) {
    dat$Q_prob[pred_q1_idx] <- predict(
      Q_reg_1,
      new_data = dat[pred_q1_idx, pi_covariates, drop = FALSE]
    )
  }

  dat$pihat <- NA_real_
  dat$pihat[Q_ind] <- dat$Q_prob[Q_ind] / (1 - dat$Q_prob[Q_ind])

  observed_idx <- dat$S == 1 & dat$R == 1 & dat$C == 0
  n_trial_hat <- c(
    sum((1 - dat$A[observed_idx]) / dat$pihat[observed_idx]),
    sum(dat$A[observed_idx] / dat$pihat[observed_idx])
  )

  n_aux <- sum(dat$wt[dat$S == 0])
  n_dat <- nrow(dat)

  ipw_0 <- ifelse(
    observed_idx & dat$A == 0,
    (dat$Y - dat$muhat_0) / dat$pihat,
    0
  )
  ipw_1 <- ifelse(
    observed_idx & dat$A == 1,
    (dat$Y - dat$muhat_1) / dat$pihat,
    0
  )
  or_0 <- ifelse(dat$S == 0, dat$muhat_0 * dat$wt, 0)
  or_1 <- ifelse(dat$S == 0, dat$muhat_1 * dat$wt, 0)

  if_0 <- (dat$S * dat$R * (1 - dat$C) * ipw_0 / n_trial_hat[1] +
             (1 - dat$S) * or_0 / n_aux) * n_dat
  if_1 <- (dat$S * dat$R * (1 - dat$C) * ipw_1 / n_trial_hat[2] +
             (1 - dat$S) * or_1 / n_aux) * n_dat

  data.frame(
    etahat_0 = mean(if_0),
    etahat_1 = mean(if_1),
    cov_00 = var(if_0) / n_dat,
    cov_01 = cov(if_0, if_1) / n_dat,
    cov_11 = var(if_1) / n_dat
  )
}


#' Prepare one reference-study replicate for HAL fitting.
#'
#' Data generation and fold construction stay in the controller process so
#' dynamic worker assignment cannot affect random-number streams. The returned
#' object contains four complete HAL job specifications plus only the analysis
#' data needed after those jobs finish.
#'
#' @param m Integer number of clusters.
#' @param n_trial Integer trial sample size.
#' @param n_aux Integer auxiliary sample size.
#' @param p_resp Numeric response probability setting.
#' @param p_cens Numeric censoring probability setting.
#' @param seed Integer simulation seed.
#' @param config List returned by `configure_fitcv_parallel()`.
#' @param cache_dir Optional cache directory for HAL fit objects.
#'
#' @return A list with named `jobs` and an `analysis` bundle.
prepare_reference_replicate <- function(m, n_trial, n_aux, p_resp, p_cens,
                                        seed, config, cache_dir = NULL) {
  sim <- generate_reference_data(
    m, n_trial, n_auxiliary = n_aux, p_resp, p_cens, seed
  )
  job_bundle <- make_reference_fit_jobs(
    dat = sim$dat,
    seed = seed,
    n_trial = sim$n_trial,
    n_aux = sim$n_aux,
    config = config,
    cache_dir = cache_dir
  )

  list(
    jobs = job_bundle$jobs,
    analysis = list(
      dat_naive = job_bundle$dat_naive,
      dat_all = job_bundle$dat_all,
      mu_covariates = job_bundle$mu_covariates,
      pi_covariates = job_bundle$pi_covariates,
      eta0 = sim$eta0,
      eta1 = sim$eta1,
      seed = as.integer(seed),
      m = as.integer(m),
      n_trial = as.integer(sim$n_trial),
      n_aux = as.integer(sim$n_aux),
      p_resp = as.numeric(p_resp),
      p_cens = as.numeric(p_cens)
    )
  )
}


#' Assemble one reference-study result after all four HAL fits finish.
#'
#' @param prepared Output from `prepare_reference_replicate()`.
#' @param fit_outputs Named outputs from `run_one_hal_job()`.
#'
#' @return A list with `results` and `fit_timings` data frames.
assemble_reference_replicate <- function(prepared, fit_outputs) {
  if (!is.list(prepared) || !is.list(prepared$analysis)) {
    stop("prepared must come from prepare_reference_replicate().", call. = FALSE)
  }
  if (!is.list(fit_outputs) || is.null(names(fit_outputs)) ||
      anyDuplicated(names(fit_outputs)) ||
      !setequal(names(fit_outputs), EXPECTED_HAL_FIT_LABELS)) {
    stop(
      "fit_outputs must contain each expected HAL fit label exactly once.",
      call. = FALSE
    )
  }
  fit_outputs <- fit_outputs[EXPECTED_HAL_FIT_LABELS]
  valid_outputs <- vapply(
    fit_outputs,
    function(x) is.list(x) && !is.null(x$fit) && is.data.frame(x$timing),
    logical(1)
  )
  if (!all(valid_outputs)) {
    stop("Every HAL fit output must contain fit and timing objects.", call. = FALSE)
  }

  analysis <- prepared$analysis

  naive_eta <- compute_naive_eta(
    dat = analysis$dat_naive,
    outcome_reg = fit_outputs$shared_outcome_reg$fit,
    censor_reg = fit_outputs$naive_censor_reg$fit,
    mu_covariates = analysis$mu_covariates,
    pi_covariates = analysis$pi_covariates
  )

  prop_eta <- compute_proposed_eta(
    dat = analysis$dat_all,
    outcome_reg = fit_outputs$shared_outcome_reg$fit,
    Q_reg_0 = fit_outputs$proposed_Q_reg_0$fit,
    Q_reg_1 = fit_outputs$proposed_Q_reg_1$fit,
    mu_covariates = analysis$mu_covariates,
    pi_covariates = analysis$pi_covariates
  )

  res <- bind_rows(
    mutate(naive_eta, name = "aipw_naive"),
    mutate(prop_eta, name = "aipw_prop")
  ) %>%
    separate(name, into = c("est", "version"), sep = "_") %>%
    mutate(
      rdhat = etahat_1 - etahat_0,
      rrhat = etahat_1 / etahat_0,
      var_rd = cov_11 + cov_00 - 2 * cov_01,
      var_rr = (cov_11 / (etahat_0^2)) +
        (cov_00 * (etahat_1^2) / (etahat_0^4)) -
        cov_01 * 2 * etahat_1 / (etahat_0^3),
      Estimator = factor(est, levels = c("aipw"), labels = c("AIPW")),
      Version = factor(
        version,
        levels = c("naive", "prop"),
        labels = c("Naive", "Proposed")
      )
    ) %>%
    select(-c(est, version)) %>%
    mutate(
      eta_0 = analysis$eta0,
      eta_1 = analysis$eta1,
      rd = eta_1 - eta_0,
      rr = eta_1 / eta_0,
      seed = analysis$seed,
      m = analysis$m,
      n_trial = analysis$n_trial,
      n_aux = analysis$n_aux,
      p_resp = analysis$p_resp,
      p_cens = analysis$p_cens
    )

  fit_timings <- bind_rows(lapply(fit_outputs, `[[`, "timing"))
  list(results = res, fit_timings = fit_timings)
}


#' Run one complete reference-study replicate with the per-run dispatcher.
#'
#' This compatibility wrapper is retained as a scientific reference path for
#' validating the global scheduler.
#'
#' @param m Integer number of clusters.
#' @param n_trial Integer trial sample size.
#' @param n_aux Integer auxiliary sample size.
#' @param p_resp Numeric response probability setting.
#' @param p_cens Numeric censoring probability setting.
#' @param seed Integer simulation seed.
#' @param config List returned by `configure_fitcv_parallel()`.
#' @param cache_dir Optional cache directory for HAL fit objects.
#' @param use_cache Logical, whether cached HAL fits may be used.
#' @param fit_cluster Optional reusable outer PSOCK cluster.
#' @param cv_cluster Optional reusable CV PSOCK cluster for serial complete
#'   fits.
#'
#' @return A list with `results` and `fit_timings` data frames.
run_reference_replicate <- function(m, n_trial, n_aux, p_resp, p_cens, seed,
                                    config, cache_dir = NULL,
                                    use_cache = TRUE, fit_cluster = NULL,
                                    cv_cluster = NULL) {
  prepared <- prepare_reference_replicate(
    m = m,
    n_trial = n_trial,
    n_aux = n_aux,
    p_resp = p_resp,
    p_cens = p_cens,
    seed = seed,
    config = config,
    cache_dir = cache_dir
  )

  fit_outputs <- run_hal_jobs(
    prepared$jobs,
    fit_workers = config$fit_workers,
    backend = config$backend,
    use_cache = use_cache,
    fit_cluster = fit_cluster,
    cv_cluster = cv_cluster
  )

  prepared$jobs <- NULL
 assemble_reference_replicate(prepared, fit_outputs)
}


#' Fit generated reference data with the formal package API.
#'
#' This adapter is the boundary between repository-only data generation and
#' the user-facing method implemented by `popart::fit_popart()`. It is useful
#' for regression comparisons against the preserved scheduler implementation.
#'
#' @param generated A list returned by `generate_reference_data()`.
#' @param control Optional object returned by `popart::popart_control()`.
#'
#' @return A `popart_fit` object.
fit_generated_data_with_popart <- function(generated, control = NULL) {
  if (!requireNamespace("popart", quietly = TRUE)) {
    stop(
      "Install the popart package before using the package-API adapter.",
      call. = FALSE
    )
  }
  if (!is.list(generated) || !is.data.frame(generated$dat) ||
      !"S" %in% names(generated$dat)) {
    stop("generated must come from generate_reference_data().", call. = FALSE)
  }

  data <- generated$dat
  covariates <- c("X1", "W1", "W2", "W3")
  trial_data <- data[data$S == 1, , drop = FALSE]
  auxiliary_data <- data[
    data$S == 0,
    c(covariates, "wt"),
    drop = FALSE
  ]
  if (is.null(control)) {
    control <- popart::popart_control()
  }

  popart::fit_popart(
    trial_data = trial_data,
    auxiliary_data = auxiliary_data,
    outcome = "Y",
    treatment = "A",
    response = "R",
    censoring = "C",
    covariates = covariates,
    auxiliary_weight = "wt",
    treatment_values = c(0, 1),
    control = control
  )
}
