###############################################################################
# PopART nuisance-fit specifications and AIPW estimating equations
###############################################################################


.POPART_CV_SEED_OFFSETS <- c(
  selection_control = 12899L,
  selection_treated = 12915L,
  outcome = 18013L,
  censoring = 14390L
)


.make_popart_nuisance_jobs <- function(data, covariates, control) {
  data$Q <- data$S * data$R * (1L - data$C)
  outcome_columns <- c("A", covariates)

  responders <- data$S == 1L & data$R == 1L
  observed <- responders & data$C == 0L
  selection_control <- (observed & data$A == 0L) | data$S == 0L
  selection_treated <- (observed & data$A == 1L) | data$S == 0L

  make_job <- function(label, X, Y, weights = NULL) {
    list(
      label = label,
      X = X,
      Y = as.numeric(Y),
      weights = if (is.null(weights)) NULL else as.numeric(weights),
      smoothness_orders = control$smoothness_orders,
      max_degree = control$max_degree,
      num_knots = control$num_knots,
      fold_ids = .make_cv_fold_ids(
        length(Y),
        control$n_cv_folds,
        control$random_seed + .POPART_CV_SEED_OFFSETS[[label]]
      ),
      n_cv_folds = control$n_cv_folds,
      n_lambda_values = control$n_lambda_values,
      n_cv_workers = control$n_cv_workers
    )
  }

  list(
    # Fixed offsets retain the validated reference fold assignments while the
    # analysis code uses domain-specific task names.
    selection_control = make_job(
      label = "selection_control",
      X = data[selection_control, covariates, drop = FALSE],
      Y = data$Q[selection_control],
      weights = data$wt[selection_control]
    ),
    selection_treated = make_job(
      label = "selection_treated",
      X = data[selection_treated, covariates, drop = FALSE],
      Y = data$Q[selection_treated],
      weights = data$wt[selection_treated]
    ),
    outcome = make_job(
      label = "outcome",
      X = data[observed, outcome_columns, drop = FALSE],
      Y = data$Y[observed]
    ),
    censoring = make_job(
      label = "censoring",
      X = data[responders, covariates, drop = FALSE],
      Y = data$C[responders]
    )
  )
}


.compute_trial_only_popart <- function(data, fits, covariates,
                                       treatment_probability) {
  data <- data[data$S == 1L & data$R == 1L, , drop = FALSE]
  outcome_columns <- c("A", covariates)

  new_control <- data[, outcome_columns, drop = FALSE]
  new_treated <- new_control
  new_control$A <- 0L
  new_treated$A <- 1L

  mu_control <- as.numeric(stats::predict(
    fits$outcome,
    new_data = new_control
  ))
  mu_treated <- as.numeric(stats::predict(
    fits$outcome,
    new_data = new_treated
  ))
  censor_probability <- as.numeric(stats::predict(
    fits$censoring,
    new_data = data[, covariates, drop = FALSE]
  ))
  censor_survival <- 1 - censor_probability
  treatment_probability_observed <- ifelse(
    data$A == 1L,
    treatment_probability,
    1 - treatment_probability
  )
  joint_probability <- treatment_probability_observed * censor_survival
  .assert_positive_predictions(
    joint_probability,
    "The trial-only treatment/censoring probability"
  )

  uncensored <- data$C == 0L
  hajek_control <- sum((1L - data$A[uncensored]) /
                         joint_probability[uncensored])
  hajek_treated <- sum(data$A[uncensored] / joint_probability[uncensored])
  .assert_hajek_denominators(hajek_control, hajek_treated, "trial-only")

  target_weight_sum <- sum(data$wt)
  n <- nrow(data)
  residual_control <- ifelse(
    uncensored & data$A == 0L,
    (data$Y - mu_control) / joint_probability,
    0
  )
  residual_treated <- ifelse(
    uncensored & data$A == 1L,
    (data$Y - mu_treated) / joint_probability,
    0
  )
  influence_control <- (
    (1L - data$C) * residual_control / hajek_control +
      mu_control * data$wt / target_weight_sum
  ) * n
  influence_treated <- (
    (1L - data$C) * residual_treated / hajek_treated +
      mu_treated * data$wt / target_weight_sum
  ) * n

  .make_popart_eta_result(
    influence_control,
    influence_treated,
    prediction_diagnostics = c(
      censor_survival_min = min(censor_survival),
      censor_survival_max = max(censor_survival),
      joint_probability_min = min(joint_probability),
      joint_probability_max = max(joint_probability)
    )
  )
}


.compute_trial_auxiliary_popart <- function(data, fits, covariates) {
  data$Q <- data$S * data$R * (1L - data$C)
  outcome_columns <- c("A", covariates)
  predict_outcome <- data$S == 0L | data$R == 1L

  new_control <- data[predict_outcome, outcome_columns, drop = FALSE]
  new_treated <- new_control
  new_control$A <- 0L
  new_treated$A <- 1L
  mu_control <- rep(NA_real_, nrow(data))
  mu_treated <- rep(NA_real_, nrow(data))
  mu_control[predict_outcome] <- as.numeric(stats::predict(
    fits$outcome,
    new_data = new_control
  ))
  mu_treated[predict_outcome] <- as.numeric(stats::predict(
    fits$outcome,
    new_data = new_treated
  ))

  observed <- data$Q == 1L
  observed_control <- observed & data$A == 0L
  observed_treated <- observed & data$A == 1L
  selection_probability <- rep(NA_real_, nrow(data))
  selection_probability[observed_control] <- as.numeric(stats::predict(
    fits$selection_control,
    new_data = data[observed_control, covariates, drop = FALSE]
  ))
  selection_probability[observed_treated] <- as.numeric(stats::predict(
    fits$selection_treated,
    new_data = data[observed_treated, covariates, drop = FALSE]
  ))
  observed_selection_probability <- selection_probability[observed]
  if (any(!is.finite(observed_selection_probability)) ||
      any(observed_selection_probability <= 0) ||
      any(observed_selection_probability >= 1)) {
    stop(
      "The estimated trial-selection probabilities must be strictly between 0 and 1.",
      call. = FALSE
    )
  }
  selection_odds <- rep(NA_real_, nrow(data))
  selection_odds[observed] <- observed_selection_probability /
    (1 - observed_selection_probability)

  hajek_control <- sum(
    (1L - data$A[observed]) / selection_odds[observed]
  )
  hajek_treated <- sum(
    data$A[observed] / selection_odds[observed]
  )
  .assert_hajek_denominators(
    hajek_control,
    hajek_treated,
    "trial-and-auxiliary"
  )

  auxiliary_weight_sum <- sum(data$wt[data$S == 0L])
  n <- nrow(data)
  residual_control <- ifelse(
    observed_control,
    (data$Y - mu_control) / selection_odds,
    0
  )
  residual_treated <- ifelse(
    observed_treated,
    (data$Y - mu_treated) / selection_odds,
    0
  )
  regression_control <- ifelse(
    data$S == 0L,
    mu_control * data$wt,
    0
  )
  regression_treated <- ifelse(
    data$S == 0L,
    mu_treated * data$wt,
    0
  )
  influence_control <- (
    data$S * data$R * (1L - data$C) * residual_control / hajek_control +
      (1L - data$S) * regression_control / auxiliary_weight_sum
  ) * n
  influence_treated <- (
    data$S * data$R * (1L - data$C) * residual_treated / hajek_treated +
      (1L - data$S) * regression_treated / auxiliary_weight_sum
  ) * n

  .make_popart_eta_result(
    influence_control,
    influence_treated,
    prediction_diagnostics = c(
      selection_probability_min = min(observed_selection_probability),
      selection_probability_max = max(observed_selection_probability),
      selection_odds_min = min(selection_odds[observed]),
      selection_odds_max = max(selection_odds[observed])
    )
  )
}


.make_popart_eta_result <- function(influence_control, influence_treated,
                                    prediction_diagnostics) {
  if (any(!is.finite(influence_control)) ||
      any(!is.finite(influence_treated))) {
    stop("Non-finite influence-function values were produced.", call. = FALSE)
  }
  n <- length(influence_control)
  covariance <- matrix(
    c(
      stats::var(influence_control) / n,
      stats::cov(influence_control, influence_treated) / n,
      stats::cov(influence_control, influence_treated) / n,
      stats::var(influence_treated) / n
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("mean_control", "mean_treated"),
      c("mean_control", "mean_treated")
    )
  )
  list(
    eta = c(
      mean_control = mean(influence_control),
      mean_treated = mean(influence_treated)
    ),
    covariance = covariance,
    influence = cbind(
      mean_control = influence_control,
      mean_treated = influence_treated
    ),
    prediction_diagnostics = prediction_diagnostics
  )
}


.assert_positive_predictions <- function(x, label) {
  if (any(!is.finite(x)) || any(x <= 0)) {
    stop(label, " contains nonpositive or non-finite values.", call. = FALSE)
  }
  invisible(x)
}


.assert_hajek_denominators <- function(control, treated, estimator) {
  if (!is.finite(control) || control <= 0 ||
      !is.finite(treated) || treated <= 0) {
    stop(
      "The ", estimator,
      " estimator has a nonpositive Hajek denominator for at least one arm.",
      call. = FALSE
    )
  }
  invisible(c(control, treated))
}


.popart_parameter_tables <- function(estimator_results, confidence_level) {
  z <- stats::qnorm(1 - (1 - confidence_level) / 2)
  rows <- lapply(names(estimator_results), function(estimator) {
    result <- estimator_results[[estimator]]
    eta <- result$eta
    covariance <- result$covariance
    estimates <- c(
      mean_control = eta[["mean_control"]],
      mean_treated = eta[["mean_treated"]],
      risk_difference = eta[["mean_treated"]] - eta[["mean_control"]],
      risk_ratio = eta[["mean_treated"]] / eta[["mean_control"]]
    )
    gradient <- rbind(
      mean_control = c(1, 0),
      mean_treated = c(0, 1),
      risk_difference = c(-1, 1),
      risk_ratio = c(
        -eta[["mean_treated"]] / eta[["mean_control"]]^2,
        1 / eta[["mean_control"]]
      )
    )
    full_covariance <- gradient %*% covariance %*% t(gradient)
    variances <- diag(full_covariance)
    standard_errors <- sqrt(pmax(variances, 0))

    list(
      estimates = data.frame(
        estimator = estimator,
        parameter = names(estimates),
        estimate = unname(estimates),
        std_error = unname(standard_errors),
        conf_low = unname(estimates - z * standard_errors),
        conf_high = unname(estimates + z * standard_errors),
        stringsAsFactors = FALSE
      ),
      variance = data.frame(
        estimator = estimator,
        parameter = names(estimates),
        variance = unname(variances),
        stringsAsFactors = FALSE
      ),
      full_covariance = full_covariance
    )
  })
  names(rows) <- names(estimator_results)

  estimates <- do.call(rbind, lapply(rows, `[[`, "estimates"))
  variance <- do.call(rbind, lapply(rows, `[[`, "variance"))
  rownames(estimates) <- NULL
  rownames(variance) <- NULL

  list(
    estimates = estimates,
    variance = variance,
    covariance = lapply(estimator_results, `[[`, "covariance"),
    full_covariance = lapply(rows, `[[`, "full_covariance")
  )
}
