#' Print a PopART fit
#'
#' Displays the four causal estimates for each estimator and version.
#'
#' @param x A `popart_fit` object returned by [fit_popart()].
#' @param digits Number of digits to display.
#' @param ... Additional arguments; currently unused.
#'
#' @return `x`, invisibly.
#'
#' @export
print.popart_fit <- function(x, digits = 3, ...) {
  estimates <- x$estimates
  output <- unique(estimates[c("estimator", "version")])
  parameters <- c(
    mean_control = "eta(0)",
    mean_treated = "eta(1)",
    risk_difference = "RD",
    risk_ratio = "RR"
  )
  output_key <- paste(output$estimator, output$version)

  for (parameter in names(parameters)) {
    rows <- estimates$parameter == parameter
    estimate_key <- paste(estimates$estimator[rows], estimates$version[rows])
    output[[parameters[[parameter]]]] <- estimates$estimate[rows][
      match(output_key, estimate_key)
    ]
  }

  cat("PopART causal estimates\n\n")
  print(output, row.names = FALSE, digits = digits)
  invisible(x)
}


#' Summarize a PopART fit
#'
#' Returns all estimator-specific causal estimates with standard errors and
#' 95 percent confidence intervals.
#'
#' @param object A `popart_fit` object returned by [fit_popart()].
#' @param digits Number of digits to display.
#' @param ... Additional arguments; currently unused.
#'
#' @return A data frame containing the estimator, version, parameter, estimate,
#'   standard error, and 95 percent confidence interval.
#'
#' @export
summary.popart_fit <- function(object, digits = 3, ...) {
  output <- object$estimates
  output[["95% CI"]] <- sprintf(
    "[%.*f, %.*f]",
    digits, output$conf_low,
    digits, output$conf_high
  )
  output$estimate <- round(output$estimate, digits)
  output$std_error <- round(output$std_error, digits)
  output <- output[c(
    "estimator", "version", "parameter", "estimate", "std_error", "95% CI"
  )]
  rownames(output) <- NULL
  output
}
