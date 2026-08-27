# Sandwich variance ---------------------------------------------------------

#' Compute a sandwich covariance matrix
#'
#' Computes an empirical sandwich covariance matrix from a stacked estimating
#' function.
#'
#' @param parameter A numeric vector containing the parameter estimates.
#' @param estimating_function A function returning observation-level estimating
#'   functions.
#' @param n A positive integer giving the sample size.
#'
#' @return The estimated covariance matrix of \code{parameter}.
#'
#' @keywords internal
#' @noRd
.sandwich_covariance <- function(parameter, estimating_function, n) {
  ## empirical bread and meat matrices
  derivative <- numDeriv::jacobian(
    function(value) colSums(estimating_function(value)),
    parameter,
    method = "simple"
  ) / -n
  psi <- estimating_function(parameter)
  meat <- crossprod(psi) / n
  inverse <- solve(derivative)
  inverse %*% meat %*% t(inverse) / n
}


#' Format potential-outcome means and covariance entries
#'
#' Combines two potential-outcome means and their covariance entries in the
#' common estimator result format.
#'
#' @param eta0 A numeric value giving the estimated mean under control.
#' @param eta1 A numeric value giving the estimated mean under treatment.
#' @param covariance A two-by-two covariance matrix for \code{eta0} and
#'   \code{eta1}.
#'
#' @return A one-row data frame containing:
#' \itemize{
#'   \item{\code{etahat_0}: estimated mean under control;}
#'   \item{\code{etahat_1}: estimated mean under treatment;}
#'   \item{\code{cov_00} and \code{cov_11}: estimated variances;}
#'   \item{\code{cov_01}: estimated covariance.}
#' }
#'
#' @keywords internal
#' @noRd
.eta_result <- function(eta0, eta1, covariance) {
  data.frame(
    etahat_0 = eta0,
    etahat_1 = eta1,
    cov_00 = covariance[1, 1],
    cov_01 = covariance[1, 2],
    cov_11 = covariance[2, 2]
  )
}
