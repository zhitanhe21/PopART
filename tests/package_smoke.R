library(popart)


expect_error <- function(expression, pattern = NULL) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  stopifnot(inherits(error, "error"))
  if (!is.null(pattern)) {
    stopifnot(grepl(pattern, conditionMessage(error), fixed = TRUE))
  }
  invisible(error)
}


exports <- getNamespaceExports("popart")
stopifnot(all(c("fit_popart", "popart_control") %in% exports))
stopifnot(!any(grepl(
  "demo|sim2|scheduler|monte.carlo|data.generation",
  exports,
  ignore.case = TRUE
)))

available_cpu_slots <- getFromNamespace(
  ".popart_available_cpu_slots",
  "popart"
)()
auto_fit_workers <- getFromNamespace(".popart_auto_fit_workers", "popart")
stopifnot(identical(
  vapply(c(1L, 2L, 3L, 4L, 5L, 6L, 8L, 64L), auto_fit_workers, integer(1)),
  c(1L, 1L, 1L, 2L, 3L, 4L, 4L, 4L)
))
expect_error(auto_fit_workers(0L), "positive integer")
expect_error(auto_fit_workers(NA_integer_), "positive integer")
expect_error(auto_fit_workers(1.5), "positive integer")

control <- popart_control()
stopifnot(
  inherits(control, "popart_control"),
  identical(control$n_cv_folds, 3L),
  identical(control$num_knots, c(5L, 3L)),
  identical(control$n_fit_workers, auto_fit_workers(available_cpu_slots)),
  identical(control$n_fit_workers_auto, TRUE),
  identical(control$detected_cpu_slots, available_cpu_slots),
  identical(control$n_cv_workers, 1L),
  identical(control$keep_nuisance_fits, FALSE)
)
auto_cv_control <- popart_control(n_cv_workers = 3L)
stopifnot(
  identical(auto_cv_control$n_fit_workers, 1L),
  identical(auto_cv_control$n_fit_workers_auto, TRUE),
  identical(auto_cv_control$n_cv_workers, 3L)
)
expect_error(
  popart_control(n_fit_workers = 2L, n_cv_workers = 2L),
  "not both"
)
expect_error(popart_control(n_cv_folds = 2L), "n_cv_folds")
expect_error(popart_control(n_fit_workers = "fast"), "auto")

base_cars <- datasets::mtcars
cars <- base_cars[
  rep(seq_len(nrow(base_cars)), 2L),
  ,
  drop = FALSE
]
rownames(cars) <- NULL
row_id <- seq_len(nrow(cars))
trial <- transform(
  cars,
  outcome = as.integer(mpg >= stats::median(mpg)),
  treatment = as.integer(row_id %% 2L == 0L),
  responded = as.integer(row_id %% 5L != 0L),
  censored = as.integer(row_id %% 3L == 0L)
)
trial$censored[trial$responded == 0L] <- NA_integer_
trial$outcome[
  trial$responded == 0L | trial$censored == 1L
] <- NA_integer_
auxiliary <- data.frame(
  wt = cars$wt,
  hp = cars$hp,
  survey_weight = 1 + (row_id %% 3L) / 10
)
stopifnot(
  nrow(trial) == 64L,
  nrow(auxiliary) == 64L,
  all(c("outcome", "treatment", "responded", "censored") %in% names(trial)),
  !any(c("outcome", "treatment", "responded", "censored") %in%
         names(auxiliary)),
  all(c("wt", "hp", "survey_weight") %in% names(auxiliary))
)

arguments <- list(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "outcome",
  treatment = "treatment",
  response = "responded",
  censoring = "censored",
  covariates = c("wt", "hp"),
  auxiliary_weight = "survey_weight",
  treatment_values = c(0, 1),
  control = control
)

prepare_input <- getFromNamespace(".prepare_popart_input", "popart")
prepared <- do.call(prepare_input, arguments)
stopifnot(
  nrow(prepared$data) == 128L,
  identical(unname(prepared$sample_sizes[c("trial", "auxiliary")]), c(64L, 64L)),
  all(prepared$data$S[seq_len(64L)] == 1L),
  all(prepared$data$S[65:128] == 0L)
)

make_jobs <- getFromNamespace(".make_popart_nuisance_jobs", "popart")
jobs <- make_jobs(prepared$data, prepared$covariates, control)
stopifnot(
  identical(
    names(jobs),
    c("selection_control", "selection_treated", "outcome", "censoring")
  ),
  all(vapply(jobs, function(job) length(job$fold_ids) == nrow(job$X), logical(1))),
  all(vapply(jobs, function(job) {
    identical(sort(unique(job$fold_ids)), seq_len(control$n_cv_folds))
  }, logical(1))),
  all(!vapply(jobs, `[[`, logical(1), "return_lasso"))
)

retained_control <- popart_control(
  n_fit_workers = 1L,
  keep_nuisance_fits = TRUE
)
retained_jobs <- make_jobs(
  prepared$data,
  prepared$covariates,
  retained_control
)
stopifnot(all(vapply(
  retained_jobs,
  `[[`,
  logical(1),
  "return_lasso"
)))

invalid_arguments <- arguments
invalid_arguments$outcome <- "missing_outcome"
expect_error(do.call(fit_popart, invalid_arguments), "missing_outcome")

if (identical(Sys.getenv("POPART_RUN_INTEGRATION_TESTS"), "true")) {
  arguments$control <- popart_control(
    n_cv_folds = 3L,
    n_lambda_values = 3L,
    num_knots = c(2L, 1L),
    n_fit_workers = 1L,
    n_cv_workers = 1L,
    keep_nuisance_fits = FALSE,
    random_seed = 20260805L
  )

  integration_input <- do.call(prepare_input, arguments)
  integration_jobs <- make_jobs(
    integration_input$data,
    integration_input$covariates,
    arguments$control
  )
  compact_job <- integration_jobs$outcome
  retained_job <- compact_job
  retained_job$return_lasso <- TRUE
  fit_hal_job <- getFromNamespace(".fit_hal_job", "popart")
  compact_hal <- fit_hal_job(compact_job)
  retained_hal <- fit_hal_job(retained_job)
  compact_prediction <- stats::predict(
    compact_hal$fit,
    new_data = compact_job$X
  )
  retained_prediction <- stats::predict(
    retained_hal$fit,
    new_data = retained_job$X
  )
  stopifnot(
    inherits(compact_hal$fit, "hal9001"),
    inherits(retained_hal$fit, "hal9001"),
    identical(names(compact_hal$fit), names(retained_hal$fit)),
    is.null(compact_hal$fit$lasso_fit),
    inherits(retained_hal$fit$lasso_fit, "cv.glmnet"),
    !is.null(compact_hal$fit$basis_list),
    !is.null(compact_hal$fit$coefs),
    !is.null(compact_hal$fit$lambda_star),
    is.matrix(compact_hal$fit$times),
    isTRUE(all.equal(
      compact_hal$fit$lambda_star,
      retained_hal$fit$lambda_star,
      tolerance = 0
    )),
    isTRUE(all.equal(
      compact_hal$fit$coefs,
      retained_hal$fit$coefs,
      tolerance = 1e-12
    )),
    isTRUE(all.equal(
      as.numeric(compact_prediction),
      as.numeric(retained_prediction),
      tolerance = 1e-12
    )),
    as.numeric(object.size(compact_hal$fit)) <
      as.numeric(object.size(retained_hal$fit))
  )

  fit <- do.call(fit_popart, arguments)
  stopifnot(
    inherits(fit, "popart_fit"),
    identical(unique(fit$estimates$estimator), c(
      "trial_only",
      "trial_auxiliary"
    )),
    identical(unique(fit$estimates$parameter), c(
      "mean_control",
      "mean_treated",
      "risk_difference",
      "risk_ratio"
    )),
    all(is.finite(fit$estimates$estimate)),
    is.null(fit$nuisance_fits),
    length(coef(fit, estimator = "trial_auxiliary")) == 4L,
    all(dim(vcov(fit, estimator = "trial_auxiliary")) == c(4L, 4L)),
    all(dim(confint(fit, estimator = "trial_auxiliary")) == c(4L, 2L))
  )
}
