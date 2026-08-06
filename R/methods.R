###############################################################################
# S3 methods for fitted PopART analyses
###############################################################################


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


#' Summarize a PopART fit
#'
#' @param object A fitted `popart_fit` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return An object of class `summary.popart_fit`.
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


#' Extract PopART coefficient estimates
#'
#' @param object A fitted `popart_fit` object.
#' @param estimator Which estimator to return: `"trial_auxiliary"`,
#'   `"trial_only"`, or `"all"`.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A named numeric vector.
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


#' Extract a PopART covariance matrix
#'
#' @param object A fitted `popart_fit` object.
#' @param estimator Which estimator covariance matrix to return.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A covariance matrix for the two arm means, risk difference, and risk
#'   ratio.
#' @export
vcov.popart_fit <- function(
    object,
    estimator = c("trial_auxiliary", "trial_only"),
    ...) {
  estimator <- match.arg(estimator)
  object$full_covariance[[estimator]]
}


#' Compute confidence intervals for a PopART fit
#'
#' @param object A fitted `popart_fit` object.
#' @param parm Optional character vector selecting parameters.
#' @param level Confidence level.
#' @param estimator Which estimator to use.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A two-column matrix of Wald confidence limits.
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


#' @export
nobs.popart_fit <- function(object, ...) {
  unname(
    object$diagnostics$sample_sizes[["trial"]] +
      object$diagnostics$sample_sizes[["auxiliary"]]
  )
}


#' @export
as.data.frame.popart_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$estimates
}
