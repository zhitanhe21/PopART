#' popart: Population-Augmented Analysis of Randomized Trials
#'
#' The `popart` package analyzes a randomized-trial sample together with an
#' auxiliary sample representing the target population. Its main interface,
#' [fit_popart()], fits trial-selection, outcome, and censoring nuisance models
#' with outer-cross-fitted highly adaptive lasso and returns augmented inverse
#' probability weighting estimates for arm-specific means, a risk difference,
#' and a risk ratio.
#'
#' @section Scope:
#' The installed package analyzes user-supplied data. Data generation, Monte
#' Carlo replication, and cross-replicate scheduling are development tools kept
#' in the repository-level `simulation/` directory and are not part of the
#' package API.
#'
#' @importFrom stats coef confint nobs vcov
#' @keywords internal
"_PACKAGE"
