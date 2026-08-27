# Parametric AIPW ----------------------------------------------------------

#' Fit a parametric AIPW estimator
#'
#' Fits either the naive or proposed parametric augmented inverse-probability
#' weighted estimator for the two mean potential outcomes.
#'
#' @param data A combined trial and auxiliary data frame containing:
#' \itemize{
#'   \item{\code{S}, \code{R}, and \code{C}: sample and observation
#'     indicators;}
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{the variables in \code{outcome_formula} and \code{weight_formula}.}
#' }
#' @param outcome_formula A logistic outcome-regression formula.
#' @param weight_formula A censoring formula for the naive estimator or combined
#'   sample-membership formula for the proposed estimator.
#' @param version A character string equal to \code{"naive"} or
#'   \code{"proposed"}.
#' @return A one-row data frame containing:
#' \itemize{
#'   \item{\code{etahat_0}: estimated mean potential outcome under control;}
#'   \item{\code{etahat_1}: estimated mean potential outcome under treatment;}
#'   \item{\code{cov_00} and \code{cov_11}: estimated variances;}
#'   \item{\code{cov_01}: estimated covariance.}
#' }
#'
#' @keywords internal
fit_aipw <- function(
    data, outcome_formula, weight_formula,
    version = c("naive", "proposed")) {
  version <- match.arg(version)
  if (version == "naive") {
    data <- data[data$S == 1L & data$R == 1L, ]
    return(.fit_aipw_naive(data, outcome_formula, weight_formula))
  }
  .fit_aipw_proposed(data, outcome_formula, weight_formula)
}
#' Fit the naive parametric AIPW estimator
#'
#' Fits the outcome and censoring regressions, constructs the naive AIPW
#' estimating functions, and computes their sandwich covariance matrix.
#'
#' @param dat A trial responder data frame containing \code{C}, \code{A},
#'   \code{Y}, \code{wt}, and the variables in \code{mu_fmla} and
#'   \code{C_fmla}.
#' @param mu_fmla A logistic outcome-regression formula.
#' @param C_fmla A logistic censoring-model formula.
#' @return A one-row data frame containing \code{etahat_0}, \code{etahat_1},
#'   \code{cov_00}, \code{cov_01}, and \code{cov_11}.
#'
#' @keywords internal
#' @noRd
.fit_aipw_naive <- function(dat, mu_fmla, C_fmla) {
  n <- nrow(dat)
  observed0 <- dat$C == 0L & dat$A == 0L
  observed1 <- dat$C == 0L & dat$A == 1L

  outcome_fit <- stats::glm(
    mu_fmla, family = stats::binomial(), data = dat[dat$C == 0L, ]
  )
  censoring_fit <- stats::glm(C_fmla, family = stats::binomial(), data = dat)

  dat0 <- dat1 <- dat
  dat0$A <- 0L
  dat1$A <- 1L
  design0 <- stats::model.matrix(mu_fmla, dat0)
  design1 <- stats::model.matrix(mu_fmla, dat1)
  design_c <- stats::model.matrix(C_fmla, dat)

  calculate <- function(beta_mu, beta_c) {
    mu0 <- stats::plogis(as.vector(design0 %*% beta_mu))
    mu1 <- stats::plogis(as.vector(design1 %*% beta_mu))
    pi <- 1 - stats::plogis(as.vector(design_c %*% beta_c))
    h0 <- sum(1 / pi[observed0])
    h1 <- sum(1 / pi[observed1])
    phi0 <- mu0
    phi1 <- mu1
    phi0[observed0] <- phi0[observed0] +
      n * (dat$Y[observed0] - mu0[observed0]) / pi[observed0] / h0
    phi1[observed1] <- phi1[observed1] +
      n * (dat$Y[observed1] - mu1[observed1]) / pi[observed1] / h1
    list(phi0 = phi0, phi1 = phi1)
  }

  beta_mu <- stats::coef(outcome_fit)
  beta_c <- stats::coef(censoring_fit)
  contribution <- calculate(beta_mu, beta_c)
  eta0 <- mean(contribution$phi0)
  eta1 <- mean(contribution$phi1)
  p_mu <- length(beta_mu)

  covariance <- .sandwich_covariance(
    c(eta0, eta1, beta_mu, beta_c),
    function(value) {
      beta_mu_value <- value[2 + seq_len(p_mu)]
      beta_c_value <- utils::tail(value, length(beta_c))
      phi <- calculate(beta_mu_value, beta_c_value)
      cbind(
        .logistic_score(dat, beta_mu_value, mu_fmla) * (1L - dat$C),
        .logistic_score(dat, beta_c_value, C_fmla),
        phi$phi0 - value[1],
        phi$phi1 - value[2]
      )
    },
    n
  )

  .eta_result(eta0, eta1, covariance)
}

#' Fit the proposed parametric AIPW estimator
#'
#' Fits the outcome and arm-specific sample-membership regressions, constructs
#' the proposed AIPW estimating functions, and computes their covariance.
#'
#' @param dat A combined trial and auxiliary data frame containing \code{S},
#'   \code{R}, \code{C}, \code{A}, \code{Y}, \code{wt}, and the variables in
#'   \code{mu_fmla} and \code{pi_fmla}.
#' @param mu_fmla A logistic outcome-regression formula.
#' @param pi_fmla A logistic sample-membership formula.
#'
#' @return A one-row data frame containing \code{etahat_0}, \code{etahat_1},
#'   \code{cov_00}, \code{cov_01}, and \code{cov_11}.
#'
#' @keywords internal
#' @noRd
.fit_aipw_proposed <- function(dat, mu_fmla, pi_fmla) {
  dat$Q <- dat$S * dat$R * (1L - dat$C)
  observed <- dat$Q == 1L
  observed0 <- observed & dat$A == 0L
  observed1 <- observed & dat$A == 1L
  auxiliary <- dat$S == 0L
  fit0 <- observed0 | auxiliary
  fit1 <- observed1 | auxiliary
  n <- nrow(dat)

  outcome_fit <- stats::glm(
    mu_fmla, family = stats::binomial(), data = dat[observed, ]
  )
  propensity0 <- do.call(stats::glm, list(
    formula = pi_fmla, family = stats::binomial(),
    data = dat[fit0, ], weights = dat$wt[fit0]
  ))
  propensity1 <- do.call(stats::glm, list(
    formula = pi_fmla, family = stats::binomial(),
    data = dat[fit1, ], weights = dat$wt[fit1]
  ))

  dat0 <- dat1 <- dat
  dat0$A <- 0L
  dat1$A <- 1L
  design0 <- stats::model.matrix(mu_fmla, dat0)
  design1 <- stats::model.matrix(mu_fmla, dat1)
  design_q <- stats::model.matrix(pi_fmla, dat)
  auxiliary_weight <- sum(dat$wt[auxiliary])

  calculate <- function(beta_mu, beta_q0, beta_q1) {
    mu0 <- stats::plogis(as.vector(design0 %*% beta_mu))
    mu1 <- stats::plogis(as.vector(design1 %*% beta_mu))
    probability0 <- stats::plogis(as.vector(design_q[observed0, , drop = FALSE] %*% beta_q0))
    probability1 <- stats::plogis(as.vector(design_q[observed1, , drop = FALSE] %*% beta_q1))
    odds0 <- probability0 / (1 - probability0)
    odds1 <- probability1 / (1 - probability1)
    h0 <- sum(1 / odds0)
    h1 <- sum(1 / odds1)

    phi0 <- phi1 <- numeric(n)
    phi0[auxiliary] <- n * dat$wt[auxiliary] * mu0[auxiliary] / auxiliary_weight
    phi1[auxiliary] <- n * dat$wt[auxiliary] * mu1[auxiliary] / auxiliary_weight
    phi0[observed0] <- n *
      (dat$Y[observed0] - mu0[observed0]) / odds0 / h0
    phi1[observed1] <- n *
      (dat$Y[observed1] - mu1[observed1]) / odds1 / h1
    list(phi0 = phi0, phi1 = phi1)
  }

  beta_mu <- stats::coef(outcome_fit)
  beta_q0 <- stats::coef(propensity0)
  beta_q1 <- stats::coef(propensity1)
  contribution <- calculate(beta_mu, beta_q0, beta_q1)
  eta0 <- mean(contribution$phi0)
  eta1 <- mean(contribution$phi1)
  p_mu <- length(beta_mu)
  p_q0 <- length(beta_q0)

  covariance <- .sandwich_covariance(
    c(eta0, eta1, beta_mu, beta_q0, beta_q1),
    function(value) {
      beta_mu_value <- value[2 + seq_len(p_mu)]
      beta_q0_value <- value[2 + p_mu + seq_len(p_q0)]
      beta_q1_value <- utils::tail(value, length(beta_q1))
      phi <- calculate(beta_mu_value, beta_q0_value, beta_q1_value)
      cbind(
        .logistic_score(dat, beta_mu_value, mu_fmla) * observed,
        .logistic_score(dat, beta_q0_value, pi_fmla) * fit0,
        .logistic_score(dat, beta_q1_value, pi_fmla) * fit1,
        phi$phi0 - value[1],
        phi$phi1 - value[2]
      )
    },
    n
  )

  .eta_result(eta0, eta1, covariance)
}
