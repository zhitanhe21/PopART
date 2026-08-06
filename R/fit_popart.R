###############################################################################
# Main analysis interface
###############################################################################


#' Fit PopART estimators to trial and auxiliary data
#'
#' Fits trial-only and trial-and-auxiliary AIPW estimators using four HAL nuisance
#' models: a shared outcome regression, one censoring regression, and separate
#' trial-selection models for the control and treated arms. This function
#' analyzes supplied data; it does not generate data or run Monte Carlo
#' simulations.
#'
#' @param trial_data A data frame containing the randomized-trial records.
#' @param auxiliary_data A data frame representing the target population. It
#'   needs the requested covariates and, optionally, an auxiliary weight column;
#'   trial-only outcome and missingness columns are not required.
#' @param outcome,treatment,response,censoring Single character strings naming
#'   columns in `trial_data`. All four variables must be binary. A value of one
#'   for `response` means the participant responded; a value of one for
#'   `censoring` means their outcome is censored.
#' @param covariates Character vector naming numeric or logical baseline
#'   covariates present in both data sets. Categorical variables should be
#'   encoded before calling this function.
#' @param auxiliary_weight Optional column in `auxiliary_data` containing
#'   nonnegative target-population weights. Unit weights are used when `NULL`.
#' @param treatment_values Two distinct values in control, treated order. The
#'   default supports a treatment column coded 0/1. For example, use
#'   `c("usual care", "intervention")` for labeled arms.
#' @param control A control object created by [popart_control()].
#'
#' @return An object of class `popart_fit` containing estimates, variances,
#'   confidence intervals, nuisance-fit timing, and diagnostic information.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_popart(
#'   trial_data = trial,
#'   auxiliary_data = auxiliary,
#'   outcome = "outcome",
#'   treatment = "arm",
#'   response = "responded",
#'   censoring = "censored",
#'   covariates = c("age", "baseline_score"),
#'   auxiliary_weight = "survey_weight"
#' )
#' summary(fit)
#' }
fit_popart <- function(
    trial_data,
    auxiliary_data,
    outcome,
    treatment,
    response,
    censoring,
    covariates,
    auxiliary_weight = NULL,
    treatment_values = c(0, 1),
    control = popart_control()) {
  if (!requireNamespace("hal9001", quietly = TRUE)) {
    stop("Package 'hal9001' is required to fit PopART models.", call. = FALSE)
  }
  control <- .validate_popart_control(control)
  input <- .prepare_popart_input(
    trial_data = trial_data,
    auxiliary_data = auxiliary_data,
    outcome = outcome,
    treatment = treatment,
    response = response,
    censoring = censoring,
    covariates = covariates,
    auxiliary_weight = auxiliary_weight,
    treatment_values = treatment_values,
    control = control
  )

  jobs <- .make_popart_nuisance_jobs(
    data = input$data,
    covariates = input$covariates,
    control = control
  )
  fit_outputs <- .fit_popart_nuisance_models(jobs, control)
  nuisance_fits <- lapply(fit_outputs, `[[`, "fit")

  estimator_results <- list(
    trial_only = .compute_trial_only_popart(
      data = input$data,
      fits = nuisance_fits,
      covariates = input$covariates,
      treatment_probability = control$treatment_probability
    ),
    trial_auxiliary = .compute_trial_auxiliary_popart(
      data = input$data,
      fits = nuisance_fits,
      covariates = input$covariates
    )
  )
  tables <- .popart_parameter_tables(
    estimator_results,
    confidence_level = control$conf_level
  )

  fit_diagnostics <- do.call(
    rbind,
    lapply(fit_outputs, function(output) output$timing)
  )
  rownames(fit_diagnostics) <- NULL
  prediction_diagnostics <- lapply(
    estimator_results,
    `[[`,
    "prediction_diagnostics"
  )
  diagnostic_messages <- .popart_diagnostic_messages(
    prediction_diagnostics,
    estimator_results,
    control$positivity_threshold
  )

  if (!isTRUE(control$keep_nuisance_fits)) {
    nuisance_fits <- NULL
  }

  structure(
    list(
      estimates = tables$estimates,
      variance = tables$variance,
      covariance = tables$covariance,
      full_covariance = tables$full_covariance,
      fit_diagnostics = fit_diagnostics,
      diagnostics = list(
        sample_sizes = input$sample_sizes,
        auxiliary_weight_scale = input$weight_scale,
        prediction_ranges = prediction_diagnostics,
        messages = diagnostic_messages
      ),
      nuisance_fits = nuisance_fits,
      columns = input$column_map,
      control = control,
      call = match.call()
    ),
    class = "popart_fit"
  )
}


.fit_popart_nuisance_models <- function(jobs, control) {
  cv_cluster <- NULL
  if (control$n_fit_workers == 1L && control$n_cv_workers > 1L) {
    cv_cluster <- .make_cv_cluster(control$n_cv_workers)
    on.exit(try(parallel::stopCluster(cv_cluster), silent = TRUE), add = TRUE)
  }

  output <- .run_nuisance_jobs(
    job_specs = jobs,
    n_fit_workers = control$n_fit_workers,
    parallel_backend = control$parallel_backend,
    cv_cluster = cv_cluster
  )
  if (!is.list(output) || !identical(names(output), names(jobs))) {
    stop("Nuisance-model fitting returned an invalid result.", call. = FALSE)
  }
  valid <- vapply(
    output,
    function(x) is.list(x) && !is.null(x$fit) && is.data.frame(x$timing),
    logical(1)
  )
  if (!all(valid)) {
    stop("At least one nuisance model did not return a valid fit.", call. = FALSE)
  }
  output
}


.popart_diagnostic_messages <- function(prediction_diagnostics,
                                        estimator_results,
                                        threshold) {
  messages <- character()
  trial_only <- prediction_diagnostics$trial_only
  trial_auxiliary <- prediction_diagnostics$trial_auxiliary

  if (trial_only[["censor_survival_min"]] < threshold) {
    messages <- c(
      messages,
      "The trial-only estimator has estimated censoring-survival probabilities near zero."
    )
  }
  if (trial_auxiliary[["selection_probability_min"]] < threshold ||
      trial_auxiliary[["selection_probability_max"]] > 1 - threshold) {
    messages <- c(
      messages,
      "The trial-and-auxiliary estimator has trial-selection probabilities near a boundary."
    )
  }
  eta_control <- vapply(
    estimator_results,
    function(x) x$eta[["mean_control"]],
    numeric(1)
  )
  if (any(abs(eta_control) < threshold)) {
    messages <- c(
      messages,
      "At least one control mean is near zero, so its risk ratio may be unstable."
    )
  }
  unique(messages)
}
