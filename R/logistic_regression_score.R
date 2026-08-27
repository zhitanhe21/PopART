# Logistic regression estimating function ---------------------------------

#' Compute logistic regression scores
#'
#' Computes the observation-level score contributions from a weighted
#' logistic regression model.
#'
#' @param data A data frame containing:
#' \itemize{
#'   \item{the response and predictors in \code{formula};}
#'   \item{\code{wt}: observation weights.}
#' }
#' @param beta A numeric vector of regression coefficients.
#' @param formula A logistic regression formula.
#'
#' @return A numeric matrix with one score row for each observation.
#'
#' @keywords internal
#' @noRd
.logistic_score <- function(data, beta, formula) {
  ## observation-level logistic regression scores
  frame <- stats::model.frame(formula, data)
  response <- stats::model.response(frame)
  design <- stats::model.matrix(formula, data)
  probability <- stats::plogis(as.vector(design %*% beta))
  design * as.vector(response - probability) * data$wt
}
