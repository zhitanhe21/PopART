# Parametric IPW ------------------------------------------------------------

#' Fit a parametric IPW estimator
#'
#' Fits either the naive or proposed parametric inverse-probability weighted
#' estimator for the two mean potential outcomes.
#'
#' @param data A combined trial and auxiliary data frame containing:
#' \itemize{
#'   \item{\code{S}: trial-membership indicator;}
#'   \item{\code{R}: response indicator;}
#'   \item{\code{C}: censoring indicator;}
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{the variables in \code{weight_formula}.}
#' }
#' @param weight_formula A censoring formula for the naive estimator or combined
#'   sample-membership formula for the proposed estimator.
#' @param version A character string equal to \code{"naive"} or
#'   \code{"proposed"}.
#' @param treatment_probability A numeric value giving the known treatment
#'   probability.
#'
#' @return A one-row data frame containing:
#' \itemize{
#'   \item{\code{etahat_0}: estimated mean potential outcome under control;}
#'   \item{\code{etahat_1}: estimated mean potential outcome under treatment;}
#'   \item{\code{cov_00} and \code{cov_11}: estimated variances;}
#'   \item{\code{cov_01}: estimated covariance.}
#' }
#'
#' @keywords internal
fit_ipw <- function(
    data, weight_formula, version = c("naive", "proposed"),
    treatment_probability = 0.5) {
  version <- match.arg(version)
  if (version == "naive") {
    data <- data[data$S == 1L & data$R == 1L, ]
    return(.fit_ipw_naive(data, weight_formula, treatment_probability))
  }
  .fit_ipw_proposed(data, weight_formula)
}
#' Fit the naive IPW calculation
#'
#' Fits the censoring model and computes the treatment-specific Hajek IPW
#' estimates and sandwich covariance matrix.
#'
#' @param data An analysis data frame containing:
#' \itemize{
#'   \item{\code{C}: censoring indicator;}
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{the variables in \code{censoring_formula}.}
#' }
#' @param censoring_formula A logistic censoring-model formula.
#' @param treatment_probability A numeric value giving the known marginal
#'   treatment probability.
#'
#' @return A one-row data frame containing:
#' \itemize{
#'   \item{\code{etahat_0}: estimated mean potential outcome under control;}
#'   \item{\code{etahat_1}: estimated mean potential outcome under treatment;}
#'   \item{\code{cov_00} and \code{cov_11}: estimated variances;}
#'   \item{\code{cov_01}: estimated covariance.}
#' }
#'
#' @keywords internal
#' @noRd
.fit_ipw_naive <- function(data, censoring_formula, treatment_probability = 0.5) {
  ## censoring and treatment probabilities
  censoring_fit <- stats::glm(
    censoring_formula, family = stats::binomial(), data = data
  )
  design <- stats::model.matrix(censoring_formula, data)
  pi_c <- 1 - stats::predict(censoring_fit, data, type = "response")
  pi_a <- ifelse(data$A == 1L, treatment_probability, 1 - treatment_probability)
  weight <- pi_c * pi_a
  observed0 <- data$C == 0L & data$A == 0L
  observed1 <- data$C == 0L & data$A == 1L

  ## Hajek IPW estimates
  denominator0 <- sum(1 / weight[observed0])
  denominator1 <- sum(1 / weight[observed1])
  eta0 <- sum(data$Y[observed0] / weight[observed0]) / denominator0
  eta1 <- sum(data$Y[observed1] / weight[observed1]) / denominator1

  ## sandwich variance estimator
  covariance <- .sandwich_covariance(
    c(eta0, eta1, denominator0, denominator1, stats::coef(censoring_fit)),
    function(value) {
      beta <- value[-seq_len(4)]
      pi_c_value <- 1 - stats::plogis(as.vector(design %*% beta))
      weight_value <- pi_c_value * pi_a
      cbind(
        .logistic_score(data, beta, censoring_formula),
        ifelse(observed0, 1 / weight_value, 0) - value[3],
        ifelse(observed1, 1 / weight_value, 0) - value[4],
        ifelse(observed0, nrow(data) * data$Y / weight_value / value[3], 0) - value[1],
        ifelse(observed1, nrow(data) * data$Y / weight_value / value[4], 0) - value[2]
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}

#' Fit the proposed IPW calculation
#'
#' Fits arm-specific sample-membership models and computes the proposed Hajek
#' IPW estimates and sandwich covariance matrix.
#'
#' @param data A combined trial and auxiliary data frame containing:
#' \itemize{
#'   \item{\code{S}, \code{R}, and \code{C}: sample and observation
#'     indicators;}
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{the variables in \code{propensity_formula}.}
#' }
#' @param propensity_formula A logistic sample-membership formula.
#'
#' @return A one-row data frame containing:
#' \itemize{
#'   \item{\code{etahat_0}: estimated mean potential outcome under control;}
#'   \item{\code{etahat_1}: estimated mean potential outcome under treatment;}
#'   \item{\code{cov_00} and \code{cov_11}: estimated variances;}
#'   \item{\code{cov_01}: estimated covariance.}
#' }
#'
#' @keywords internal
#' @noRd
.fit_ipw_proposed <- function(data, propensity_formula) {
  ## observed-outcome and auxiliary analysis sample
  data$Q <- data$S * data$R * (1L - data$C)
  restricted <- data$Q == 1L | data$S == 0L
  fit0 <- restricted & (data$A == 0L | data$Q == 0L)
  fit1 <- restricted & (data$A == 1L | data$Q == 0L)
  data0 <- data[fit0, ]
  data1 <- data[fit1, ]

  ## treatment-specific sample-membership models
  propensity0 <- do.call(stats::glm, list(
    formula = propensity_formula,
    family = stats::binomial(),
    data = data0,
    weights = data0$wt
  ))
  propensity1 <- do.call(stats::glm, list(
    formula = propensity_formula,
    family = stats::binomial(),
    data = data1,
    weights = data1$wt
  ))
  design <- stats::model.matrix(propensity_formula, data)
  observed0 <- data$Q == 1L & data$A == 0L
  observed1 <- data$Q == 1L & data$A == 1L

  ## estimated selection odds
  probability <- rep(0, nrow(data))
  probability[observed0] <- stats::predict(propensity0, data[observed0, ], type = "response")
  probability[observed1] <- stats::predict(propensity1, data[observed1, ], type = "response")
  odds <- probability / (1 - probability)

  ## Hajek IPW estimates
  denominator0 <- sum(1 / odds[observed0])
  denominator1 <- sum(1 / odds[observed1])
  eta0 <- sum(data$Y[observed0] / odds[observed0]) / denominator0
  eta1 <- sum(data$Y[observed1] / odds[observed1]) / denominator1
  length0 <- length(stats::coef(propensity0))

  ## sandwich variance estimator
  covariance <- .sandwich_covariance(
    c(
      eta0, eta1, denominator0, denominator1,
      stats::coef(propensity0), stats::coef(propensity1)
    ),
    function(value) {
      beta0 <- value[4 + seq_len(length0)]
      beta1 <- utils::tail(value, length(stats::coef(propensity1)))
      probability_value <- rep(0, nrow(data))
      probability_value[observed0] <- stats::plogis(
        as.vector(design[observed0, , drop = FALSE] %*% beta0)
      )
      probability_value[observed1] <- stats::plogis(
        as.vector(design[observed1, , drop = FALSE] %*% beta1)
      )
      odds_value <- probability_value / (1 - probability_value)
      cbind(
        .logistic_score(data, beta0, propensity_formula) * fit0,
        .logistic_score(data, beta1, propensity_formula) * fit1,
        ifelse(observed0, 1 / odds_value, 0) - value[3],
        ifelse(observed1, 1 / odds_value, 0) - value[4],
        ifelse(observed0, nrow(data) * data$Y / odds_value / value[3], 0) - value[1],
        ifelse(observed1, nrow(data) * data$Y / odds_value / value[4], 0) - value[2]
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}
