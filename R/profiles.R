###############################################################################
# Convenience run profiles for the local Simulation 2 workflow
###############################################################################


SIM2_PROFILE_NAMES <- c("smoke", "quick", "hour", "formal")


#' Normalize and validate a Simulation 2 profile name.
#'
#' @param profile Character profile name.
#'
#' @return One lower-case profile name.
normalize_sim2_profile <- function(profile) {
  if (length(profile) != 1L || is.na(profile)) {
    stop(
      "profile must be one of: ",
      paste(SIM2_PROFILE_NAMES, collapse = ", "),
      call. = FALSE
    )
  }

  profile <- tolower(trimws(as.character(profile)))
  if (!nzchar(profile) || !profile %in% SIM2_PROFILE_NAMES) {
    stop(
      "profile must be one of: ",
      paste(SIM2_PROFILE_NAMES, collapse = ", "),
      call. = FALSE
    )
  }
  profile
}


#' Return the convenience defaults for one Simulation 2 profile.
#'
#' Profiles are argument presets only. They are intentionally not included in
#' fit-cache identities or run-checkpoint signatures.
#'
#' @param profile Character profile name.
#'
#' @return A named list of arguments accepted by `run_sim2_demo()`.
#' @export
sim2_profile_defaults <- function(profile) {
  profile <- normalize_sim2_profile(profile)

  common_small <- list(
    scale = "small",
    max_runs = 4L,
    reserve_cores = 0L,
    fit_workers = 2L,
    cv_workers_per_fit = 1L,
    cv_folds = SIM2_DEFAULT_CV_FOLDS,
    cv_nlambda = 50L,
    backend = "psock",
    use_cache = FALSE
  )

  switch(
    profile,
    smoke = list(
      scale = "demo",
      max_runs = 1L,
      mc_reps = 1L,
      reserve_cores = 0L,
      fit_workers = 1L,
      cv_workers_per_fit = 1L,
      cv_folds = SIM2_DEFAULT_CV_FOLDS,
      cv_nlambda = 50L,
      backend = "auto",
      use_cache = FALSE
    ),
    quick = utils::modifyList(
      common_small,
      list(mc_reps = 1L),
      keep.null = TRUE
    ),
    hour = utils::modifyList(
      common_small,
      list(mc_reps = 2L),
      keep.null = TRUE
    ),
    formal = utils::modifyList(
      common_small,
      list(mc_reps = 20L),
      keep.null = TRUE
    )
  )
}


#' Overlay explicit arguments on a Simulation 2 profile.
#'
#' Resolution order is legacy defaults, then profile defaults, then explicit
#' overrides. `NULL` is retained so API callers can explicitly override a
#' profile value such as `max_runs`.
#'
#' @param profile Character profile name.
#' @param overrides Named list of explicit user arguments.
#' @param legacy_defaults Optional named list of caller-specific legacy
#'   defaults.
#'
#' @return A named list of resolved arguments.
resolve_sim2_profile_args <- function(profile, overrides = list(),
                                      legacy_defaults = list()) {
  if (!is.list(overrides) || !is.list(legacy_defaults)) {
    stop("overrides and legacy_defaults must be lists.", call. = FALSE)
  }
  if (length(overrides) > 0L &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("Every profile override must be named.", call. = FALSE)
  }
  if (length(legacy_defaults) > 0L &&
      (is.null(names(legacy_defaults)) ||
       any(!nzchar(names(legacy_defaults))))) {
    stop("Every legacy default must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(overrides))) {
    stop("Profile overrides must not contain duplicate names.", call. = FALSE)
  }

  resolved <- utils::modifyList(
    legacy_defaults,
    sim2_profile_defaults(profile),
    keep.null = TRUE
  )
  utils::modifyList(resolved, overrides, keep.null = TRUE)
}


#' Read one profile default while preserving a caller-specific legacy default.
#'
#' @param defaults Named list returned by `sim2_profile_defaults()`, or an
#'   empty list when no profile was requested.
#' @param name Character argument name.
#' @param legacy Value used when the profile does not define `name`.
#'
#' @return The profile value when present, otherwise `legacy`.
sim2_profile_default <- function(defaults, name, legacy) {
  if (name %in% names(defaults)) defaults[[name]] else legacy
}


format_sim2_profile_value <- function(value) {
  if (is.null(value)) return("NULL")
  if (length(value) == 0L) return("<empty>")
  paste(as.character(value), collapse = ",")
}


#' Print the resolved arguments supplied by a convenience profile.
#'
#' Effective worker counts are printed later by `run_sim2_demo()` after local
#' CPU-budget clipping.
#' @noRd
print_sim2_profile_resolution <- function(profile, resolved_args,
                                          explicit_names = character()) {
  profile <- normalize_sim2_profile(profile)
  if (!is.list(resolved_args)) {
    stop("resolved_args must be a list.", call. = FALSE)
  }

  cat("Simulation 2 convenience profile\n")
  cat("  profile = ", profile, "\n", sep = "")
  for (name in names(resolved_args)) {
    source_label <- if (name %in% explicit_names) " [explicit]" else ""
    cat(
      "  ",
      name,
      " = ",
      format_sim2_profile_value(resolved_args[[name]]),
      source_label,
      "\n",
      sep = ""
    )
  }
  cat(
    "  effective CPU/worker values are shown in Local parallel plan\n"
  )
  invisible(NULL)
}


#' Run Simulation 2 through a named convenience profile.
#'
#' The wrapper leaves `run_sim2_demo()` and its defaults unchanged. Named
#' arguments in `...` override the selected profile.
#'
#' @param project_dir Optional source-checkout path. Installed package users can
#'   leave it as `NULL`.
#' @param profile Character profile name: `smoke`, `quick`, `hour`, or
#'   `formal`.
#' @param ... Named arguments forwarded to `run_sim2_demo()`.
#'
#' @return Invisibly returns the output from `run_sim2_demo()` with profile
#'   metadata added to the returned list.
#' @export
run_sim2_demo_profile <- function(project_dir = NULL, profile, ...) {
  if (!exists("run_sim2_demo", mode = "function", inherits = TRUE)) {
    stop(
      "run_sim2_demo() is not available in the current session.",
      call. = FALSE
    )
  }

  profile <- normalize_sim2_profile(profile)
  overrides <- list(...)
  if (length(overrides) > 0L &&
      (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("All run_sim2_demo_profile() overrides must be named.", call. = FALSE)
  }
  if (anyDuplicated(names(overrides))) {
    stop("Profile overrides must not contain duplicate names.", call. = FALSE)
  }

  allowed <- setdiff(names(formals(run_sim2_demo)), "project_dir")
  unknown <- setdiff(names(overrides), allowed)
  if (length(unknown) > 0L) {
    stop(
      "Unknown run_sim2_demo() override(s): ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  resolved <- resolve_sim2_profile_args(
    profile,
    overrides = overrides
  )
  print_sim2_profile_resolution(
    profile,
    resolved,
    explicit_names = names(overrides)
  )

  call_args <- c(list(project_dir = project_dir), resolved)
  output <- do.call(run_sim2_demo, call_args)
  output$profile <- profile
  output$profile_arguments <- resolved
  invisible(output)
}
