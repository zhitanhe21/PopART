# Parametric G-formula ------------------------------------------------------

#' Fit a parametric G-formula estimator
#'
#' Fits either the naive or proposed parametric G-formula estimator for the two
#' mean potential outcomes.
#'
#' @param data A combined trial and auxiliary data frame containing:
#' \itemize{
#'   \item{\code{S}: trial-membership indicator;}
#'   \item{\code{R}: response indicator;}
#'   \item{\code{C}: censoring indicator;}
#'   \item{\code{A}: treatment indicator;}
#'   \item{\code{Y}: outcome;}
#'   \item{\code{wt}: analysis weight;}
#'   \item{the variables in \code{outcome_formula}.}
#' }
#' @param outcome_formula A logistic outcome-regression formula.
#' @param version A character string equal to \code{"naive"} or
#'   \code{"proposed"}.
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
fit_gformula <- function(data, outcome_formula, version = c("naive", "proposed")) {
  version <- match.arg(version)
  proposed <- version == "proposed"

  if (!proposed) {
    data <- data[data$S == 1L & data$R == 1L, ]
  }

  fit_rows <- if (proposed) {
    data$S == 1L & data$R == 1L & data$C == 0L
  } else {
    data$C == 0L
  }
  target_rows <- if (proposed) data$S == 0L else rep(TRUE, nrow(data))

  outcome_fit <- stats::glm(
    outcome_formula, family = stats::binomial(), data = data[fit_rows, ]
  )

  data0 <- data1 <- data
  data0$A <- 0L
  data1$A <- 1L
  design0 <- stats::model.matrix(outcome_formula, data0)
  design1 <- stats::model.matrix(outcome_formula, data1)
  mu0 <- stats::predict(outcome_fit, data0, type = "response")
  mu1 <- stats::predict(outcome_fit, data1, type = "response")
  eta0 <- stats::weighted.mean(mu0[target_rows], data$wt[target_rows])
  eta1 <- stats::weighted.mean(mu1[target_rows], data$wt[target_rows])

  covariance <- .sandwich_covariance(
    c(eta0, eta1, stats::coef(outcome_fit)),
    function(value) {
      beta <- value[-c(1, 2)]
      cbind(
        .logistic_score(data, beta, outcome_formula) * as.numeric(fit_rows),
        as.numeric(target_rows) * data$wt *
          (stats::plogis(as.vector(design0 %*% beta)) - value[1]),
        as.numeric(target_rows) * data$wt *
          (stats::plogis(as.vector(design1 %*% beta)) - value[2])
      )
    },
    nrow(data)
  )

  .eta_result(eta0, eta1, covariance)
}
