###############################################################################
# Main analysis interface
###############################################################################


#' Fit population-augmented randomized-trial estimators
#'
#' Fits augmented inverse probability weighting (AIPW) estimators to a
#' randomized-trial sample and an auxiliary sample from the target population.
#' Both a trial-only estimator and a trial-and-auxiliary estimator are returned.
#'
#' Four highly adaptive lasso (HAL) nuisance models are fitted: one outcome
#' regression, one censoring regression, and arm-specific trial-selection
#' regressions. The outcome regression is shared by the two estimators. The
#' function analyzes data supplied by the user; it does not generate data or
#' run Monte Carlo simulations.
#'
#' The trial-and-auxiliary AIPW framework follows Richardson, Shook-Sa, and
#' Hudgens (n.d.). The current use of HAL for nuisance fitting is an
#' implementation choice of this package and is cited separately below.
#'
#' `trial_data` and `auxiliary_data` must contain the same baseline covariates,
#' but the auxiliary sample does not need treatment, response, censoring, or
#' outcome columns. The outcome may be missing only when it is unobserved
#' (`response = 0` or `censoring = 1`). Censoring may be missing when
#' `response = 0`; such entries are treated as zero internally. Input rows are
#' not modified.
#'
#' The returned Wald intervals use the influence-function covariance estimates
#' and the confidence level in `control`. They are not truncated to the natural
#' parameter space. Likewise, estimated nuisance probabilities are diagnosed
#' but are not truncated.
#'
#' @param trial_data A `data.frame` or tibble with one row per randomized-trial
#'   participant. It must contain the four trial variables named below and every
#'   column in `covariates`.
#' @param auxiliary_data A `data.frame` or tibble with one row per auxiliary
#'   target-population record. It must contain every column in `covariates` and,
#'   when requested, the column named by `auxiliary_weight`.
#' @param outcome A nonempty character string naming a binary 0/1 or
#'   `FALSE`/`TRUE` outcome column in `trial_data`.
#' @param treatment A nonempty character string naming the treatment column in
#'   `trial_data`. Its two allowed values are supplied in control-then-treated
#'   order through `treatment_values`.
#' @param response A nonempty character string naming a binary 0/1 or
#'   `FALSE`/`TRUE` response-indicator column in `trial_data`; one means that the
#'   participant responded.
#' @param censoring A nonempty character string naming a binary 0/1 or
#'   `FALSE`/`TRUE` censoring-indicator column in `trial_data`; one means that a
#'   responding participant's outcome is censored.
#' @param covariates A nonempty character vector of distinct baseline-covariate
#'   column names present in both data sets. Covariates must be numeric or
#'   logical and cannot contain missing or non-finite values. Encode categorical
#'   variables before fitting.
#' @param auxiliary_weight `NULL` or a nonempty character string naming a column
#'   of finite, nonnegative target-population weights in `auxiliary_data`. The
#'   weights must have a positive sum. Unit weights are used when `NULL`.
#' @param treatment_values A length-two vector containing the treatment values
#'   in control-then-treated order. The default, `c(0, 1)`, supports numeric 0/1
#'   coding; labeled arms may instead use, for example,
#'   `c("usual care", "intervention")`.
#' @param control A `popart_control` object created by [popart_control()].
#'
#' @return An object of class `popart_fit`, implemented as a list with the
#'   following elements:
#'   \itemize{
#'   \item `estimates`: a data frame with columns `estimator`, `parameter`,
#'     `estimate`, `std_error`, `conf_low`, and `conf_high`.
#'   \item `variance`: a data frame with columns `estimator`, `parameter`, and
#'     `variance` for each arm mean, risk difference, and risk ratio.
#'   \item `covariance`: a named list of two-by-two covariance matrices for the
#'     control and treated arm means, one matrix per estimator.
#'   \item `full_covariance`: a named list of four-by-four covariance matrices
#'     for the arm means, risk difference, and risk ratio, one per estimator.
#'   \item `fit_diagnostics`: a data frame with columns `fit`, `observations`,
#'     `predictors`, `response_mean`, `weight_sum`, `elapsed_seconds`,
#'     `basis_seconds`, `design_matrix_seconds`, `lasso_seconds`,
#'     `selected_lambda`, `n_cv_folds`, `n_lambda_values`, `n_cv_workers`,
#'     `process_id`, `started_at`, and `finished_at`.
#'   \item `diagnostics`: a list containing `sample_sizes`, the auxiliary-weight
#'     normalization scale (`auxiliary_weight_scale`), nuisance-prediction ranges
#'     (`prediction_ranges`), and diagnostic `messages`.
#'   \item `nuisance_fits`: a list named `selection_control`,
#'     `selection_treated`, `outcome`, and `censoring` containing fitted HAL
#'     objects when `control$keep_nuisance_fits` is `TRUE`; otherwise `NULL`.
#'   \item `columns`: a list with elements `outcome`, `treatment`, `response`,
#'     `censoring`, `covariates`, `auxiliary_weight`, and `treatment_values`,
#'     recording the input-column specification used for the analysis.
#'   \item `control`: the validated `popart_control` object used for fitting.
#'   \item `call`: the matched function call.
#'   }
#'
#' @references
#' Richardson, B. D., Shook-Sa, B. E., and Hudgens, M. G. (n.d.).
#' "Causal Inference from Cluster-Randomized Trials with Differential
#' Nonresponse." Unpublished manuscript submitted to *Biometrics*.
#'
#' Hejazi, N. S., Coyle, J. R., and van der Laan, M. J. (2020).
#' "hal9001: Scalable highly adaptive lasso regression in R."
#' *Journal of Open Source Software*. \doi{10.21105/joss.02526}.
#'
#' @examples
#' trial_file <- system.file("extdata", "example_trial.csv", package = "popart")
#' auxiliary_file <- system.file(
#'   "extdata", "example_auxiliary.csv", package = "popart"
#' )
#' trial <- utils::read.csv(trial_file)
#' auxiliary <- utils::read.csv(auxiliary_file)
#'
#' # Reduced HAL settings keep the help-page example quick.
#' example_control <- popart_control(
#'   n_lambda_values = 5L,
#'   max_degree = 1L,
#'   num_knots = 1L
#' )
#' fit <- fit_popart(
#'   trial_data = trial,
#'   auxiliary_data = auxiliary,
#'   outcome = "outcome",
#'   treatment = "treatment",
#'   response = "responded",
#'   censoring = "censored",
#'   covariates = c("baseline_risk", "baseline_binary"),
#'   auxiliary_weight = "survey_weight",
#'   control = example_control
#' )
#' fit$estimates
#'
#' @export
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
