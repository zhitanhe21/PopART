###############################################################################
# S3 methods for fitted PopART analyses
###############################################################################


#' Methods for fitted PopART analyses
#'
#' Standard S3 methods extract, summarize, print, and reformat objects returned
#' by [fit_popart()]. The print method shows the trial-and-auxiliary estimates;
#' [summary()] and [as.data.frame()] retain results for both the trial-only and
#' trial-and-auxiliary estimators.
#'
#' `coef()` returns point estimates, `vcov()` returns their joint covariance
#' matrix, and `confint()` computes normal-approximation Wald intervals. The
#' covariance matrix is ordered as `mean_control`, `mean_treated`,
#' `risk_difference`, and `risk_ratio`. `nobs()` counts both trial and auxiliary
#' records because both samples contribute to the population-augmented analysis.
#'
#' @param x A `popart_fit` object for `print()` and `as.data.frame()`, or a
#'   `summary.popart_fit` object for the summary print method.
#' @param object A `popart_fit` object returned by [fit_popart()].
#' @param digits A single integer controlling the number of significant digits
#'   printed. The default is derived from `getOption("digits")`.
#' @param estimator A character string selecting an estimator. `coef()` accepts
#'   `"trial_auxiliary"`, `"trial_only"`, or `"all"`; `vcov()` and `confint()`
#'   accept `"trial_auxiliary"` or `"trial_only"`.
#' @param parm `NULL` or a character vector selecting any of `"mean_control"`,
#'   `"mean_treated"`, `"risk_difference"`, and `"risk_ratio"`. `NULL` selects
#'   all four parameters.
#' @param level A single numeric confidence level strictly between zero and one.
#'   By default, the level used to fit `object` is used.
#' @param row.names `NULL` or a value supplied for compatibility with
#'   `as.data.frame()`; currently ignored.
#' @param optional A logical value supplied for compatibility with
#'   `as.data.frame()`; currently ignored.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The methods return the following objects:
#'   \itemize{
#'   \item `print.popart_fit()` and `print.summary.popart_fit()` return their
#'     input object invisibly after displaying it.
#'   \item `summary.popart_fit()` returns a `summary.popart_fit` list containing
#'     `call`, `estimates`, `sample_sizes`, `fit_diagnostics`, diagnostic
#'     `messages`, and `conf_level`.
#'   \item `coef.popart_fit()` returns a named numeric vector of estimates.
#'   \item `vcov.popart_fit()` returns a four-by-four numeric covariance matrix.
#'   \item `confint.popart_fit()` returns a numeric matrix with one row per
#'     selected parameter and two columns containing the lower and upper Wald
#'     confidence limits.
#'   \item `nobs.popart_fit()` returns the total number of trial and auxiliary
#'     records as an unnamed scalar.
#'   \item `as.data.frame.popart_fit()` returns the `estimates` data frame with
#'     estimator, parameter, estimate, standard-error, and interval columns.
#'   }
#'
#' @examples
#' trial <- utils::read.csv(system.file(
#'   "extdata", "example_trial.csv", package = "popart"
#' ))
#' auxiliary <- utils::read.csv(system.file(
#'   "extdata", "example_auxiliary.csv", package = "popart"
#' ))
#' fit <- fit_popart(
#'   trial_data = trial,
#'   auxiliary_data = auxiliary,
#'   outcome = "outcome",
#'   treatment = "treatment",
#'   response = "responded",
#'   censoring = "censored",
#'   covariates = c("baseline_risk", "baseline_binary"),
#'   auxiliary_weight = "survey_weight",
#'   control = popart_control(
#'     n_lambda_values = 5L,
#'     max_degree = 1L,
#'     num_knots = 1L
#'   )
#' )
#'
#' print(fit)
#' summary(fit)
#' coef(fit, estimator = "trial_auxiliary")
#' vcov(fit, estimator = "trial_auxiliary")
#' confint(fit, parm = c("risk_difference", "risk_ratio"))
#' nobs(fit)
#' head(as.data.frame(fit))
#'
#' @name popart_fit_methods
NULL


#' @rdname popart_fit_methods
#' @export
print.popart_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("PopART fit\n")
  cat(
    "  Trial / auxiliary observations: ",
    x$diagnostics$sample_sizes[["trial"]], " / ",
    x$diagnostics$sample_sizes[["auxiliary"]], "\n",
    sep = ""
  )
  displayed <- x$estimates[
    x$estimates$estimator == "trial_auxiliary",
    c("parameter", "estimate", "std_error", "conf_low", "conf_high"),
    drop = FALSE
  ]
  print(displayed, row.names = FALSE, digits = digits)
  if (length(x$diagnostics$messages)) {
    cat("  Diagnostic message(s): ", length(x$diagnostics$messages), "\n", sep = "")
  }
  invisible(x)
}


#' @rdname popart_fit_methods
#' @export
summary.popart_fit <- function(object, ...) {
  structure(
    list(
      call = object$call,
      estimates = object$estimates,
      sample_sizes = object$diagnostics$sample_sizes,
      fit_diagnostics = object$fit_diagnostics,
      messages = object$diagnostics$messages,
      conf_level = object$control$conf_level
    ),
    class = "summary.popart_fit"
  )
}


#' @rdname popart_fit_methods
#' @export
print.summary.popart_fit <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...) {
  cat("PopART analysis summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nSample sizes:\n")
  print(x$sample_sizes)
  cat("\nEstimates and ", 100 * x$conf_level, "% Wald intervals:\n", sep = "")
  print(x$estimates, row.names = FALSE, digits = digits)
  if (length(x$messages)) {
    cat("\nDiagnostics:\n")
    for (message in x$messages) {
      cat("- ", message, "\n", sep = "")
    }
  }
  invisible(x)
}


#' @rdname popart_fit_methods
#' @export
coef.popart_fit <- function(
    object,
    estimator = c("trial_auxiliary", "trial_only", "all"),
    ...) {
  estimator <- match.arg(estimator)
  table <- object$estimates
  if (!identical(estimator, "all")) {
    table <- table[table$estimator == estimator, , drop = FALSE]
  }
  names <- if (identical(estimator, "all")) {
    paste(table$estimator, table$parameter, sep = ":")
  } else {
    table$parameter
  }
  stats::setNames(table$estimate, names)
}


#' @rdname popart_fit_methods
#' @export
vcov.popart_fit <- function(
    object,
    estimator = c("trial_auxiliary", "trial_only"),
    ...) {
  estimator <- match.arg(estimator)
  object$full_covariance[[estimator]]
}


#' @rdname popart_fit_methods
#' @export
confint.popart_fit <- function(
    object,
    parm = NULL,
    level = object$control$conf_level,
    estimator = c("trial_auxiliary", "trial_only"),
    ...) {
  estimator <- match.arg(estimator)
  level <- .popart_probability(level, "level")
  table <- object$estimates[
    object$estimates$estimator == estimator,
    c("parameter", "estimate"),
    drop = FALSE
  ]
  covariance <- object$full_covariance[[estimator]]
  if (!is.null(parm)) {
    if (!is.character(parm) || anyNA(parm)) {
      stop("parm must be a character vector of parameter names.", call. = FALSE)
    }
    unknown <- setdiff(parm, table$parameter)
    if (length(unknown)) {
      stop("Unknown parameter(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    table <- table[match(parm, table$parameter), , drop = FALSE]
    covariance <- covariance[parm, parm, drop = FALSE]
  }
  standard_error <- sqrt(pmax(diag(covariance), 0))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  result <- cbind(
    lower = table$estimate - critical * standard_error,
    upper = table$estimate + critical * standard_error
  )
  rownames(result) <- table$parameter
  colnames(result) <- c(
    paste0(format(100 * (1 - level) / 2, trim = TRUE), "%"),
    paste0(format(100 * (1 + level) / 2, trim = TRUE), "%")
  )
  result
}


#' @rdname popart_fit_methods
#' @export
nobs.popart_fit <- function(object, ...) {
  unname(
    object$diagnostics$sample_sizes[["trial"]] +
      object$diagnostics$sample_sizes[["auxiliary"]]
  )
}


#' @rdname popart_fit_methods
#' @export
as.data.frame.popart_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$estimates
}
