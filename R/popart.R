#' popart: Population-Augmented Analysis of Randomized Trials
#'
#' The package fits naive and proposed g-formula, IPW, and AIPW estimators.
#'
#' @keywords internal
"_PACKAGE"


#' Fit naive and proposed PopART estimators
#'
#' Fits naive and proposed G-formula, IPW, and AIPW estimators by combining a
#' randomized trial with an auxiliary target-population sample. Four shared
#' HAL nuisance models are fitted once and reused across the six estimators.
#'
#' The trial and auxiliary samples are converted to a common internal data
#' structure before the nuisance models and estimators are fitted.
#'
#' @param trial_data A data frame containing the trial observations and the
#'   following variables:
#' \itemize{
#'   \item{the column named by \code{outcome}: binary outcome;}
#'   \item{the column named by \code{treatment}: treatment assignment;}
#'   \item{the column named by \code{response}: binary response indicator;}
#'   \item{the column named by \code{censoring}: binary censoring indicator;}
#'   \item{the columns named by \code{covariates}: baseline covariates.}
#' }
#' @param auxiliary_data A data frame containing:
#' \itemize{
#'   \item{the columns named by \code{covariates}: baseline covariates;}
#'   \item{optionally, the column named by \code{auxiliary_weight}: sampling
#'     weights.}
#' }
#' @param outcome A character string naming the binary outcome column.
#' @param treatment A character string naming the treatment column.
#' @param response A character string naming the binary response column.
#' @param censoring A character string naming the binary censoring column, where
#'   1 denotes a censored outcome.
#' @param cluster A character string naming the cluster identifier column in
#'   both samples. HAL folds and estimated variances are grouped by cluster.
#' @param covariates A character vector naming the baseline covariates used by
#'   the HAL nuisance models.
#' @param auxiliary_weight An optional character string naming the auxiliary
#'   sampling-weight column. If \code{NULL}, all auxiliary observations receive
#'   weight one.
#' @param treatment_values A vector containing the control and treated values,
#'   in that order.
#' @param treatment_probability A numeric value giving the known marginal
#'   probability of assignment to the treated value.
#' @param n_cv_folds An integer giving the number of HAL cross-validation folds.
#' @param n_lambda_values An integer giving the number of lasso penalties.
#' @param random_seed An integer seed used to construct the HAL folds.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item{\code{estimates}: estimates, standard errors, and 95 percent
#'     confidence intervals for both potential-outcome means, the risk
#'     difference, and the risk ratio;}
#'   \item{\code{variance}: the corresponding estimated variances;}
#'   \item{\code{covariance}: six named covariance matrices for the two
#'     estimated potential-outcome means.}
#' }
#'
#' @examples
#' \dontrun{
#' trial <- data.frame(
#'   Y = c(0, 1, 1, 0), A = c(0, 0, 1, 1),
#'   R = 1, C = 0, cluster = c(1, 1, 2, 2), X = c(0, 0, 1, 1)
#' )
#' auxiliary <- data.frame(
#'   cluster = c(1, 1, 2, 2), X = c(0, 0, 1, 1)
#' )
#'
#' fit <- fit_popart(
#'   trial, auxiliary,
#'   outcome = "Y", treatment = "A", response = "R", censoring = "C",
#'   cluster = "cluster", covariates = "X"
#' )
#' fit$estimates
#' }
#'
#' @export
fit_popart <- function(
    trial_data,
    auxiliary_data,
    outcome,
    treatment,
    response,
    censoring,
    cluster,
    covariates,
    auxiliary_weight = NULL,
    treatment_values = c(0, 1),
    treatment_probability = 0.5,
    n_cv_folds = 5L,
    n_lambda_values = 50L,
    random_seed = 1L) {

  # Prepare the combined trial and auxiliary data.
  data <- .prepare_popart_data(
    trial_data, auxiliary_data, outcome, treatment, response, censoring,
    cluster, covariates, auxiliary_weight, treatment_values
  )
  jobs <- .make_popart_jobs(
    data, covariates, n_cv_folds, n_lambda_values, random_seed, data$.cluster
  )

  # Fit the four shared HAL nuisance models.
  fits <- lapply(jobs, .fit_hal_job)

  # Compute all six estimators from the shared fits.
  fit <- .finish_popart(
    data, covariates, fits, treatment_probability, data$.cluster
  )
  class(fit) <- "popart_fit"
  fit
}


#' Prepare the combined PopART analysis data
#'
#' Renames the analysis variables, normalizes auxiliary sampling weights, and
#' combines the trial and auxiliary samples in the internal data format.
#'
#' @param trial_data A trial data frame containing:
#' \itemize{
#'   \item{the columns named by \code{outcome}, \code{treatment},
#'     \code{response}, and \code{censoring};}
#'   \item{the columns named by \code{covariates}.}
#' }
#' @param auxiliary_data An auxiliary data frame containing:
#' \itemize{
#'   \item{the columns named by \code{covariates};}
#'   \item{when supplied, the column named by \code{auxiliary_weight}.}
#' }
#' @param outcome A character string naming the outcome column.
#' @param treatment A character string naming the treatment column.
#' @param response A character string naming the response column.
#' @param censoring A character string naming the censoring column.
#' @param cluster A character string naming the cluster identifier column.
#' @param covariates A character vector naming the baseline covariates.
#' @param auxiliary_weight An optional character string naming the auxiliary
#'   weight column.
#' @param treatment_values A vector containing the control and treated values.
#'
#' @return A combined data frame containing:
#' \itemize{
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{R}: response indicator;}
#'   \item{\code{C}: censoring indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{\code{S}: trial-membership indicator;}
#'   \item{\code{.cluster}: cluster identifier;}
#'   \item{the baseline covariates.}
#' }
#'
#' @keywords internal
#' @noRd
.prepare_popart_data <- function(
    trial_data, auxiliary_data, outcome, treatment, response, censoring,
    cluster, covariates, auxiliary_weight = NULL,
    treatment_values = c(0, 1)) {
  A <- match(as.character(trial_data[[treatment]]), treatment_values) - 1L
  R <- as.integer(trial_data[[response]])
  C <- as.integer(trial_data[[censoring]])
  Y <- as.numeric(trial_data[[outcome]])
  C[is.na(C) & R == 0L] <- 0L
  Y[is.na(Y) & (R == 0L | C == 1L)] <- 0
  trial_x <- trial_data[covariates]
  auxiliary_x <- auxiliary_data[covariates]
  trial_cluster <- trial_data[[cluster]]
  auxiliary_cluster <- auxiliary_data[[cluster]]

  auxiliary_weights <- if (is.null(auxiliary_weight)) {
    rep(1, nrow(auxiliary_data))
  } else {
    as.numeric(auxiliary_data[[auxiliary_weight]])
  }
  auxiliary_weights <- auxiliary_weights / mean(auxiliary_weights)

  trial <- cbind(
    data.frame(
      A = A, R = R, C = C, Y = Y, wt = 1, S = 1L,
      .cluster = trial_cluster
    ),
    trial_x
  )
  auxiliary <- cbind(
    data.frame(
      A = 0L, R = 0L, C = 0L, Y = 0,
      wt = auxiliary_weights, S = 0L, .cluster = auxiliary_cluster
    ),
    auxiliary_x
  )
  rbind(trial, auxiliary)
}


#' Create the four shared HAL fitting jobs
#'
#' Creates the arm-specific selection jobs, outcome-regression job, and
#' censoring-model job used by all six estimators.
#'
#' @param data A combined analysis data frame returned by
#'   \code{.prepare_popart_data()}.
#' @param covariates A character vector naming the baseline covariates.
#' @param folds Number of cross-validation folds.
#' @param lambdas Number of lasso penalties.
#' @param seed Random-number seed used to construct the folds.
#' @param cluster Cluster identifier for each row of \code{data}.
#'
#' @return A named list containing four HAL fitting jobs.
#'
#' @keywords internal
#' @noRd
.make_popart_jobs <- function(
    data, covariates, folds, lambdas, seed,
    cluster = seq_len(nrow(data))) {
  data$Q <- data$S * data$R * (1L - data$C)
  responders <- data$S == 1L & data$R == 1L
  observed <- responders & data$C == 0L
  q0 <- (observed & data$A == 0L) | data$S == 0L
  q1 <- (observed & data$A == 1L) | data$S == 0L

  job <- function(rows, predictors, outcome, weights = NULL, offset) {
    set.seed(seed + offset)
    groups <- unique(cluster[rows])
    group_fold <- sample(rep(seq_len(folds), length.out = length(groups)))
    list(
      X = data[rows, predictors, drop = FALSE],
      Y = as.numeric(outcome[rows]),
      weights = if (is.null(weights)) NULL else weights[rows],
      foldid = group_fold[match(cluster[rows], groups)],
      folds = folds,
      lambdas = lambdas
    )
  }

  list(
    selection_control = job(q0, covariates, data$Q, data$wt, 12899L),
    selection_treated = job(q1, covariates, data$Q, data$wt, 12915L),
    outcome = job(observed, c("A", covariates), data$Y, offset = 18013L),
    censoring = job(responders, covariates, data$C, offset = 14390L)
  )
}


#' Fit one HAL nuisance model
#'
#' Fits one job created by \code{.make_popart_jobs()} without starting another
#' parallel backend.
#'
#' @param job A list containing the predictors, outcome, optional weights,
#'   cross-validation folds, and lambda settings.
#'
#' @return A fitted model from \code{hal9001::fit_hal()}.
#'
#' @keywords internal
#' @noRd
.fit_hal_job <- function(job) {
  hal9001::fit_hal(
    X = job$X,
    Y = job$Y,
    weights = job$weights,
    family = "binomial",
    smoothness_orders = 1,
    max_degree = 2,
    num_knots = c(5L, 3L),
    return_lasso = FALSE,
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-4,
      nfolds = job$folds,
      nlambda = job$lambdas,
      foldid = job$foldid,
      parallel = FALSE,
      prediction_bounds = "default"
    )
  )
}


#' Compute the six estimators from fitted HAL models
#'
#' Uses the four shared nuisance fits to compute naive and proposed G-formula,
#' IPW, and AIPW contributions.
#'
#' @param data A combined analysis data frame.
#' @param covariates A character vector naming the baseline covariates.
#' @param fits A named list containing the four fitted HAL models.
#' @param treatment_probability Known probability of treatment assignment.
#' @param cluster Cluster identifier for each row of \code{data}.
#'
#' @return Estimate, variance, and covariance tables for all six estimators.
#'
#' @keywords internal
#' @noRd
.finish_popart <- function(
    data, covariates, fits, treatment_probability = 0.5,
    cluster = seq_len(nrow(data))) {
  n <- nrow(data)
  responders <- data$S == 1L & data$R == 1L
  auxiliary <- data$S == 0L
  observed <- responders & data$C == 0L
  observed0 <- observed & data$A == 0L
  observed1 <- observed & data$A == 1L

  # Outcome predictions ---------------------------------------------------

  predict_rows <- responders | auxiliary
  new0 <- new1 <- data[predict_rows, c("A", covariates), drop = FALSE]
  new0$A <- 0L
  new1$A <- 1L
  mu0 <- mu1 <- numeric(n)
  mu0[predict_rows] <- stats::predict(fits$outcome, new_data = new0)
  mu1[predict_rows] <- stats::predict(fits$outcome, new_data = new1)

  # Naive and proposed weights -------------------------------------------

  pi_c <- rep(NA_real_, n)
  pi_c[responders] <- 1 - stats::predict(
    fits$censoring, new_data = data[responders, covariates, drop = FALSE]
  )
  pi_naive <- ifelse(
    data$A == 1L, treatment_probability, 1 - treatment_probability
  ) * pi_c

  selection_probability <- rep(NA_real_, n)
  selection_probability[observed0] <- stats::predict(
    fits$selection_control,
    new_data = data[observed0, covariates, drop = FALSE]
  )
  selection_probability[observed1] <- stats::predict(
    fits$selection_treated,
    new_data = data[observed1, covariates, drop = FALSE]
  )
  odds <- selection_probability / (1 - selection_probability)

  n_response <- sum(responders)
  auxiliary_weight <- sum(data$wt[auxiliary])
  h_naive_0 <- sum(1 / pi_naive[observed0])
  h_naive_1 <- sum(1 / pi_naive[observed1])
  h_proposed_0 <- sum(1 / odds[observed0])
  h_proposed_1 <- sum(1 / odds[observed1])

  # Estimating-function contributions ------------------------------------

  contributions <- list(
    gformula_naive = cbind(
      n * ifelse(responders, mu0 / n_response, 0),
      n * ifelse(responders, mu1 / n_response, 0)
    ),
    gformula_proposed = cbind(
      n * ifelse(auxiliary, data$wt * mu0 / auxiliary_weight, 0),
      n * ifelse(auxiliary, data$wt * mu1 / auxiliary_weight, 0)
    ),
    ipw_naive = cbind(
      n * ifelse(observed0, data$Y / pi_naive / h_naive_0, 0),
      n * ifelse(observed1, data$Y / pi_naive / h_naive_1, 0)
    ),
    ipw_proposed = cbind(
      n * ifelse(observed0, data$Y / odds / h_proposed_0, 0),
      n * ifelse(observed1, data$Y / odds / h_proposed_1, 0)
    ),
    aipw_naive = cbind(
      n * (
        ifelse(observed0, (data$Y - mu0) / pi_naive / h_naive_0, 0) +
          ifelse(responders, mu0 / n_response, 0)
      ),
      n * (
        ifelse(observed1, (data$Y - mu1) / pi_naive / h_naive_1, 0) +
          ifelse(responders, mu1 / n_response, 0)
      )
    ),
    aipw_proposed = cbind(
      n * (
        ifelse(observed0, (data$Y - mu0) / odds / h_proposed_0, 0) +
          ifelse(auxiliary, data$wt * mu0 / auxiliary_weight, 0)
      ),
      n * (
        ifelse(observed1, (data$Y - mu1) / odds / h_proposed_1, 0) +
          ifelse(auxiliary, data$wt * mu1 / auxiliary_weight, 0)
      )
    )
  )

  .format_popart_results(contributions, cluster)
}


#' Format PopART estimator results
#'
#' Converts the six pairs of arm-specific contributions into estimate,
#' variance, and covariance tables.
#'
#' @param contributions A named list of two-column contribution matrices.
#' @param cluster Cluster identifier for each contribution row. Rows from the
#'   same cluster are summed before the covariance is calculated.
#'
#' @return A list containing \code{estimates}, \code{variance}, and
#'   \code{covariance}.
#'
#' @keywords internal
#' @noRd
.format_popart_results <- function(
    contributions,
    cluster = seq_len(nrow(contributions[[1]]))) {
  output <- lapply(names(contributions), function(name) {
    contribution <- contributions[[name]]
    eta <- colMeans(contribution)
    centered <- sweep(contribution, 2L, eta)
    cluster_sum <- rowsum(centered, cluster, reorder = FALSE)
    n <- nrow(contribution)
    n_clusters <- nrow(cluster_sum)
    covariance <- crossprod(cluster_sum) * n_clusters /
      ((n_clusters - 1) * n^2)
    dimnames(covariance) <- list(
      c("mean_control", "mean_treated"),
      c("mean_control", "mean_treated")
    )

    estimate <- c(
      mean_control = eta[[1]],
      mean_treated = eta[[2]],
      risk_difference = eta[[2]] - eta[[1]],
      risk_ratio = eta[[2]] / eta[[1]]
    )
    gradient <- rbind(
      mean_control = c(1, 0),
      mean_treated = c(0, 1),
      risk_difference = c(-1, 1),
      risk_ratio = c(-eta[[2]] / eta[[1]]^2, 1 / eta[[1]])
    )
    variance <- diag(gradient %*% covariance %*% t(gradient))
    standard_error <- sqrt(variance)
    label <- strsplit(name, "_", fixed = TRUE)[[1]]
    estimator <- c(
      gformula = "G-Formula", ipw = "IPW", aipw = "AIPW"
    )[[label[[1]]]]
    version <- c(naive = "Naive", proposed = "Proposed")[[label[[2]]]]

    list(
      estimates = data.frame(
        estimator = estimator,
        version = version,
        parameter = names(estimate),
        estimate = unname(estimate),
        std_error = unname(standard_error),
        conf_low = unname(estimate - 1.96 * standard_error),
        conf_high = unname(estimate + 1.96 * standard_error)
      ),
      variance = data.frame(
        estimator = estimator,
        version = version,
        parameter = names(estimate),
        variance = unname(variance)
      ),
      covariance = covariance
    )
  })

  list(
    estimates = do.call(rbind, lapply(output, `[[`, "estimates")),
    variance = do.call(rbind, lapply(output, `[[`, "variance")),
    covariance = stats::setNames(
      lapply(output, `[[`, "covariance"), names(contributions)
    )
  )
}
