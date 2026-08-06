###############################################################################
# Input validation and canonical analysis data
###############################################################################


.prepare_popart_input <- function(
    trial_data,
    auxiliary_data,
    outcome,
    treatment,
    response,
    censoring,
    covariates,
    auxiliary_weight,
    treatment_values,
    control) {
  .assert_popart_data_frame(trial_data, "trial_data")
  .assert_popart_data_frame(auxiliary_data, "auxiliary_data")

  roles <- c(
    outcome = .assert_popart_column(outcome, "outcome"),
    treatment = .assert_popart_column(treatment, "treatment"),
    response = .assert_popart_column(response, "response"),
    censoring = .assert_popart_column(censoring, "censoring")
  )
  covariates <- .assert_popart_covariates(covariates)
  auxiliary_weight <- .assert_popart_column(
    auxiliary_weight,
    "auxiliary_weight",
    allow_null = TRUE
  )

  named_roles <- c(unname(roles), covariates)
  if (anyDuplicated(named_roles)) {
    stop(
      "outcome, treatment, response, censoring, and covariates must refer to ",
      "distinct columns.",
      call. = FALSE
    )
  }
  if (!is.null(auxiliary_weight) && auxiliary_weight %in% covariates) {
    stop("auxiliary_weight cannot also be a covariate.", call. = FALSE)
  }

  .require_popart_columns(trial_data, named_roles, "trial_data")
  auxiliary_columns <- c(covariates, auxiliary_weight)
  .require_popart_columns(auxiliary_data, auxiliary_columns, "auxiliary_data")

  treatment_values <- .validate_treatment_values(treatment_values)
  trial_treatment <- .encode_treatment(
    trial_data[[roles[["treatment"]]]],
    treatment_values,
    roles[["treatment"]]
  )
  trial_response <- .encode_binary(
    trial_data[[roles[["response"]]]],
    roles[["response"]],
    allow_missing = FALSE
  )
  trial_censoring <- .encode_binary(
    trial_data[[roles[["censoring"]]]],
    roles[["censoring"]],
    allow_missing = TRUE
  )

  missing_censoring <- is.na(trial_censoring)
  if (any(missing_censoring & trial_response == 1L)) {
    stop(
      "censoring cannot be missing among trial participants with response = 1.",
      call. = FALSE
    )
  }
  trial_censoring[missing_censoring & trial_response == 0L] <- 0L
  if (any(trial_response == 0L & trial_censoring == 1L)) {
    stop(
      "censoring must be 0 or missing when response = 0.",
      call. = FALSE
    )
  }

  trial_outcome <- .encode_binary(
    trial_data[[roles[["outcome"]]]],
    roles[["outcome"]],
    allow_missing = TRUE
  )
  observed <- trial_response == 1L & trial_censoring == 0L
  if (anyNA(trial_outcome[observed])) {
    stop(
      "outcome cannot be missing where response = 1 and censoring = 0.",
      call. = FALSE
    )
  }
  trial_outcome[!observed & is.na(trial_outcome)] <- 0L

  x_names <- sprintf("x%03d", seq_along(covariates))
  trial_x <- .numeric_covariate_frame(
    trial_data,
    covariates,
    x_names,
    "trial_data"
  )
  auxiliary_x <- .numeric_covariate_frame(
    auxiliary_data,
    covariates,
    x_names,
    "auxiliary_data"
  )

  if (is.null(auxiliary_weight)) {
    raw_auxiliary_weight <- rep(1, nrow(auxiliary_data))
  } else {
    raw_auxiliary_weight <- suppressWarnings(
      as.numeric(auxiliary_data[[auxiliary_weight]])
    )
    if (length(raw_auxiliary_weight) != nrow(auxiliary_data) ||
        any(!is.finite(raw_auxiliary_weight)) ||
        any(raw_auxiliary_weight < 0) ||
        sum(raw_auxiliary_weight) <= 0) {
      stop(
        "auxiliary_weight must contain finite nonnegative values with a ",
        "positive sum.",
        call. = FALSE
      )
    }
  }
  auxiliary_weight_scale <- 1
  normalized_auxiliary_weight <- raw_auxiliary_weight
  if (isTRUE(control$normalize_auxiliary_weights)) {
    auxiliary_weight_scale <- mean(raw_auxiliary_weight)
    if (!is.finite(auxiliary_weight_scale) || auxiliary_weight_scale <= 0) {
      stop("The mean auxiliary weight must be positive.", call. = FALSE)
    }
    normalized_auxiliary_weight <- raw_auxiliary_weight / auxiliary_weight_scale
  }

  trial_internal <- data.frame(
    A = trial_treatment,
    R = trial_response,
    C = trial_censoring,
    Y = trial_outcome,
    wt = 1,
    S = 1L,
    trial_x,
    check.names = FALSE
  )
  auxiliary_internal <- data.frame(
    A = 0L,
    R = 0L,
    C = 0L,
    Y = 0L,
    wt = normalized_auxiliary_weight,
    S = 0L,
    auxiliary_x,
    check.names = FALSE
  )
  combined <- rbind(trial_internal, auxiliary_internal)
  rownames(combined) <- NULL

  .validate_popart_fit_samples(combined, x_names, control$n_cv_folds)

  list(
    data = combined,
    covariates = x_names,
    column_map = list(
      outcome = unname(roles[["outcome"]]),
      treatment = unname(roles[["treatment"]]),
      response = unname(roles[["response"]]),
      censoring = unname(roles[["censoring"]]),
      covariates = stats::setNames(x_names, covariates),
      auxiliary_weight = auxiliary_weight,
      treatment_values = treatment_values
    ),
    weight_scale = auxiliary_weight_scale,
    sample_sizes = c(
      trial = nrow(trial_data),
      auxiliary = nrow(auxiliary_data),
      trial_responders = sum(trial_response == 1L),
      trial_observed_outcomes = sum(observed),
      observed_control = sum(observed & trial_treatment == 0L),
      observed_treated = sum(observed & trial_treatment == 1L)
    )
  )
}


.assert_popart_data_frame <- function(x, name) {
  if (!is.data.frame(x)) {
    stop(name, " must be a data.frame or tibble.", call. = FALSE)
  }
  if (nrow(x) < 1L) {
    stop(name, " must contain at least one row.", call. = FALSE)
  }
  if (is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) ||
      anyDuplicated(names(x))) {
    stop(name, " must have unique, nonempty column names.", call. = FALSE)
  }
  invisible(x)
}


.assert_popart_column <- function(x, argument, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(argument, " must be one nonempty column name.", call. = FALSE)
  }
  x
}


.assert_popart_covariates <- function(x) {
  if (!is.character(x) || length(x) < 1L || anyNA(x) ||
      any(!nzchar(x)) || anyDuplicated(x)) {
    stop(
      "covariates must be a nonempty character vector of unique column names.",
      call. = FALSE
    )
  }
  x
}


.require_popart_columns <- function(data, columns, data_name) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) {
    stop(
      data_name, " is missing required column", if (length(missing_columns) > 1L) "s" else "",
      ": ", paste(missing_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(data)
}


.validate_treatment_values <- function(x) {
  if (length(x) != 2L || anyNA(x)) {
    stop(
      "treatment_values must contain exactly two values in control, treated order.",
      call. = FALSE
    )
  }
  value <- as.character(x)
  if (anyDuplicated(value)) {
    stop("treatment_values must contain two distinct values.", call. = FALSE)
  }
  stats::setNames(value, c("control", "treated"))
}


.encode_treatment <- function(x, treatment_values, column) {
  if (anyNA(x)) {
    stop("Treatment column '", column, "' cannot contain missing values.", call. = FALSE)
  }
  encoded <- match(as.character(x), unname(treatment_values)) - 1L
  if (anyNA(encoded)) {
    unexpected <- unique(as.character(x)[is.na(encoded)])
    stop(
      "Treatment column '", column, "' contains values outside treatment_values: ",
      paste(utils::head(unexpected, 5L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!all(c(0L, 1L) %in% encoded)) {
    stop("Both treatment groups must be represented in trial_data.", call. = FALSE)
  }
  encoded
}


.encode_binary <- function(x, column, allow_missing) {
  if (is.logical(x)) {
    value <- as.integer(x)
  } else {
    character_value <- as.character(x)
    value <- rep(NA_integer_, length(character_value))
    value[character_value == "0"] <- 0L
    value[character_value == "1"] <- 1L
  }
  invalid <- is.na(value) & !is.na(x)
  if (any(invalid)) {
    stop("Column '", column, "' must be coded as 0/1 or FALSE/TRUE.", call. = FALSE)
  }
  if (!isTRUE(allow_missing) && anyNA(value)) {
    stop("Column '", column, "' cannot contain missing values.", call. = FALSE)
  }
  value
}


.numeric_covariate_frame <- function(data, columns, internal_names, data_name) {
  out <- lapply(columns, function(column) {
    value <- data[[column]]
    if (!is.numeric(value) && !is.logical(value)) {
      stop(
        "Covariate '", column, "' in ", data_name,
        " must be numeric or logical; encode categorical variables before fitting.",
        call. = FALSE
      )
    }
    value <- as.numeric(value)
    if (any(!is.finite(value))) {
      stop(
        "Covariate '", column, "' in ", data_name,
        " cannot contain missing or non-finite values.",
        call. = FALSE
      )
    }
    value
  })
  names(out) <- internal_names
  as.data.frame(out, check.names = FALSE)
}


.validate_popart_fit_samples <- function(data, covariates, n_cv_folds) {
  responders <- data$S == 1L & data$R == 1L
  observed <- responders & data$C == 0L

  samples <- list(
    outcome = list(index = observed, response = data$Y[observed]),
    censoring = list(index = responders, response = data$C[responders]),
    selection_control = list(
      index = (observed & data$A == 0L) | data$S == 0L,
      response = data$S[(observed & data$A == 0L) | data$S == 0L]
    ),
    selection_treated = list(
      index = (observed & data$A == 1L) | data$S == 0L,
      response = data$S[(observed & data$A == 1L) | data$S == 0L]
    )
  )

  for (name in names(samples)) {
    sample <- samples[[name]]
    n <- sum(sample$index)
    if (n < n_cv_folds) {
      stop(
        "The ", name, " nuisance fit has ", n,
        " observations, fewer than n_cv_folds = ", n_cv_folds, ".",
        call. = FALSE
      )
    }
    if (length(unique(sample$response)) < 2L) {
      stop(
        "The ", name, " nuisance response must contain both 0 and 1.",
        call. = FALSE
      )
    }
  }
  invisible(covariates)
}
