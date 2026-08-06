###############################################################################
# Safe hashing, cache, checkpoint, and atomic-output helpers
###############################################################################


FIT_CACHE_SCHEMA <- "popart_reference_fit_cache_v1"
RUN_CHECKPOINT_SCHEMA <- "popart_reference_run_checkpoint_v1"
REFERENCE_STUDY_IMPLEMENTATION <- "popart_reference_study_v1"
REFERENCE_STUDY_CODE_SCHEMA <- "popart_reference_algorithm_20260805_v1"
EXPECTED_HAL_FIT_LABELS <- c(
  "proposed_Q_reg_0",
  "proposed_Q_reg_1",
  "shared_outcome_reg",
  "naive_censor_reg"
)


#' Return a package version as a plain character value.
#'
#' @param package Character package name.
#'
#' @return Character package version, or `NA_character_` when unavailable.
package_version_string <- function(package) {
  tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
}


#' Hash an R object without requiring a hard dependency on `digest`.
#'
#' Uses SHA-256 through `digest` when available. The base-R fallback serializes
#' to a temporary RDS file and calculates an MD5 checksum.
#'
#' @param object Any serializable R object.
#'
#' @return Character content hash.
hash_r_object <- function(object) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(object, algo = "sha256", serialize = TRUE))
  }

  temp_file <- tempfile(pattern = "aipw_hash_", fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)
  saveRDS(object, temp_file, compress = FALSE, version = 3)
  unname(tools::md5sum(temp_file))
}


#' Make a string safe for use as one file-name component.
#' @noRd
safe_path_component <- function(value) {
  gsub("[^A-Za-z0-9_.-]", "_", as.character(value))
}


#' Create a fingerprint for the computational reference-study source code.
#'
#' @param project_dir Character path to the repository root.
#'
#' @return Character hash based on the R files that can affect simulation
#'   results or checkpoint contents. Tests, analysis scripts, and reporting
#'   helpers are intentionally excluded so presentation-only edits do not
#'   invalidate otherwise compatible checkpoints.
make_code_fingerprint <- function(project_dir = NULL) {
  if (is.null(project_dir) || !nzchar(project_dir)) {
    project_dir <- system.file(package = "popart")
  }
  if (!nzchar(project_dir)) {
    return(hash_r_object(list(
      schema = REFERENCE_STUDY_CODE_SCHEMA,
      package_version = package_version_string("popart")
    )))
  }

  project_dir <- normalizePath(project_dir)
  relative <- c(
    file.path("simulation", "R", "configuration.R"),
    file.path("simulation", "R", "persistence.R"),
    file.path("simulation", "R", "data_generation.R"),
    file.path("simulation", "R", "hal_fitting.R"),
    file.path("simulation", "R", "hal_scheduler.R"),
    file.path("simulation", "R", "reference_estimators.R"),
    file.path("simulation", "R", "monte_carlo.R")
  )
  files <- file.path(project_dir, relative)
  keep <- file.exists(files)
  files <- files[keep]
  relative <- relative[keep]
  if (length(files) == 0L) {
    return(hash_r_object(list(
      schema = REFERENCE_STUDY_CODE_SCHEMA,
      package_version = package_version_string("popart")
    )))
  }

  checksums <- unname(tools::md5sum(files))
  hash_r_object(list(
    schema = REFERENCE_STUDY_CODE_SCHEMA,
    source_checksums = stats::setNames(checksums, relative)
  ))
}


#' Build the complete statistical identity for one HAL fit.
#'
#' Resource-layout values such as backend and worker counts are deliberately
#' excluded. They change scheduling, not the fitted statistical task.
#' @noRd
make_fit_cache_identity <- function(label, X, Y, weights, family,
                                    smoothness_orders, max_degree, num_knots,
                                    foldid, cv_folds, cv_nlambda) {
  manifest <- list(
    schema = FIT_CACHE_SCHEMA,
    label = label,
    X = X,
    Y = Y,
    weights = weights,
    family = family,
    smoothness_orders = smoothness_orders,
    max_degree = max_degree,
    num_knots = num_knots,
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-4,
      nfolds = cv_folds,
      nlambda = cv_nlambda,
      foldid = foldid,
      prediction_bounds = "default"
    ),
    package_versions = list(
      hal9001 = package_version_string("hal9001"),
      glmnet = package_version_string("glmnet"),
      Matrix = package_version_string("Matrix")
    ),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform
  )

  list(
    schema = FIT_CACHE_SCHEMA,
    key = hash_r_object(manifest),
    package_versions = manifest$package_versions,
    R_version = manifest$R_version,
    platform = manifest$platform
  )
}


#' Build a v2 cache path for one HAL fit.
#' @noRd
make_fit_cache_file <- function(cache_dir, label, seed, n_trial, n_aux,
                                cache_key) {
  if (is.null(cache_dir) || !nzchar(cache_dir) ||
      is.null(cache_key) || !nzchar(cache_key)) {
    return(NULL)
  }

  label_safe <- safe_path_component(label)
  label_dir <- file.path(cache_dir, "v2", label_safe)
  dir.create(label_dir, recursive = TRUE, showWarnings = FALSE)
  file_name <- paste0(
    "seed", seed,
    "_ntrial", n_trial,
    "_naux", n_aux,
    "_", cache_key,
    ".rds"
  )
  file.path(label_dir, file_name)
}


#' Move an invalid file aside without deleting it.
#' @noRd
quarantine_invalid_file <- function(path, reason = "invalid") {
  if (is.null(path) || !file.exists(path)) {
    return(invisible(NULL))
  }

  suffix <- paste0(
    ".", safe_path_component(reason), "-",
    format(Sys.time(), "%Y%m%d%H%M%S"), "-pid", Sys.getpid()
  )
  quarantine_path <- paste0(path, suffix)
  moved <- file.rename(path, quarantine_path)
  if (isTRUE(moved)) {
    warning(
      "Moved an invalid file aside: ", quarantine_path,
      call. = FALSE
    )
    return(invisible(quarantine_path))
  }

  warning(
    "Could not quarantine invalid file; it will be ignored: ", path,
    call. = FALSE
  )
  invisible(NULL)
}


#' Promote a same-directory temporary file to its final path.
#' @noRd
recover_atomic_target <- function(target_file) {
  if (file.exists(target_file)) {
    return(invisible(FALSE))
  }

  backups <- Sys.glob(paste0(target_file, ".previous-*"))
  if (length(backups) == 0L) {
    return(invisible(FALSE))
  }
  info <- file.info(backups)
  backup_file <- backups[[which.max(info$mtime)]]
  recovered <- file.rename(backup_file, target_file)
  if (!isTRUE(recovered)) {
    warning(
      "Could not recover prior atomic-output backup: ", backup_file,
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  warning(
    "Recovered interrupted atomic-output replacement: ", target_file,
    call. = FALSE
  )
  invisible(TRUE)
}


#' Promote a same-directory temporary file to its final path.
#' @noRd
atomic_promote_file <- function(temp_file, target_file, overwrite = FALSE) {
  if (!file.exists(temp_file)) {
    stop("Temporary file does not exist: ", temp_file, call. = FALSE)
  }
  dir.create(dirname(target_file), recursive = TRUE, showWarnings = FALSE)
  recover_atomic_target(target_file)

  if (!file.exists(target_file)) {
    promoted <- file.rename(temp_file, target_file)
    if (!isTRUE(promoted) && !overwrite && file.exists(target_file)) {
      return(invisible(FALSE))
    }
    if (!isTRUE(promoted)) {
      stop("Could not promote temporary file to: ", target_file, call. = FALSE)
    }
    return(invisible(TRUE))
  }

  if (!isTRUE(overwrite)) {
    return(invisible(FALSE))
  }

  backup_file <- paste0(
    target_file,
    ".previous-", format(Sys.time(), "%Y%m%d%H%M%S"),
    "-pid", Sys.getpid()
  )
  if (!file.rename(target_file, backup_file)) {
    stop("Could not prepare existing file for replacement: ", target_file,
         call. = FALSE)
  }

  promoted <- file.rename(temp_file, target_file)
  if (!isTRUE(promoted)) {
    restored <- file.rename(backup_file, target_file)
    if (!isTRUE(restored)) {
      stop(
        "Could not replace or restore target. Backup remains at: ",
        backup_file,
        call. = FALSE
      )
    }
    stop("Could not atomically replace: ", target_file, call. = FALSE)
  }

  unlink(backup_file)
  invisible(TRUE)
}


#' Save an R object through a same-directory temporary file.
#' @noRd
atomic_save_rds <- function(object, path, compress = FALSE,
                            overwrite = FALSE, verify = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !isTRUE(overwrite)) {
    return(invisible(FALSE))
  }

  temp_file <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".partial"
  )
  on.exit(unlink(temp_file), add = TRUE)
  saveRDS(object, temp_file, compress = compress, version = 3)

  if (is.na(file.info(temp_file)$size) || file.info(temp_file)$size < 1L) {
    stop("Temporary RDS was not written correctly for: ", path, call. = FALSE)
  }
  if (isTRUE(verify)) {
    verified <- tryCatch({
      readRDS(temp_file)
      TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(verified)) {
      stop("Could not verify temporary RDS file for: ", path, call. = FALSE)
    }
  }

  atomic_promote_file(temp_file, path, overwrite = overwrite)
}


#' Write a data frame through a same-directory temporary CSV.
#' @noRd
atomic_write_csv <- function(data, path, overwrite = TRUE) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path) && !isTRUE(overwrite)) {
    return(invisible(FALSE))
  }

  temp_file <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".partial"
  )
  on.exit(unlink(temp_file), add = TRUE)
  utils::write.csv(data, temp_file, row.names = FALSE)

  if (is.na(file.info(temp_file)$size) || file.info(temp_file)$size < 1L) {
    stop("Temporary CSV was not written correctly for: ", path, call. = FALSE)
  }
  verified <- tryCatch({
    utils::read.csv(temp_file, check.names = FALSE)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(verified)) {
    stop("Could not verify temporary CSV file for: ", path, call. = FALSE)
  }
  atomic_promote_file(temp_file, path, overwrite = overwrite)
}


#' Validate and read one fit-cache object.
#' @noRd
read_fit_cache <- function(path, expected_key) {
  if (is.null(path)) {
    return(NULL)
  }
  recover_atomic_target(path)
  if (!file.exists(path)) {
    return(NULL)
  }

  cached <- tryCatch(readRDS(path), error = function(e) e)
  valid <- !inherits(cached, "error") &&
    is.list(cached) &&
    !is.null(cached$fit) &&
    is.data.frame(cached$timing) &&
    is.list(cached$cache_meta) &&
    identical(cached$cache_meta$schema, FIT_CACHE_SCHEMA) &&
    identical(cached$cache_meta$key, expected_key)

  if (!isTRUE(valid)) {
    quarantine_invalid_file(path, "invalid-fit-cache")
    return(NULL)
  }
  cached
}


#' Create the signature for one complete simulation run.
#' @noRd
make_run_signature <- function(project_dir, seed, m, n_trial, n_aux,
                               p_resp, p_cens, config,
                               code_fingerprint = NULL) {
  if (is.null(code_fingerprint)) {
    code_fingerprint <- make_code_fingerprint(project_dir)
  }
  manifest <- list(
    schema = RUN_CHECKPOINT_SCHEMA,
    implementation = REFERENCE_STUDY_IMPLEMENTATION,
    code_fingerprint = code_fingerprint,
    seed = as.integer(seed),
    m = as.integer(m),
    requested_n_trial = as.integer(n_trial),
    requested_n_aux = as.integer(n_aux),
    p_resp = as.numeric(p_resp),
    p_cens = as.numeric(p_cens),
    hal_controls = list(
      family = "binomial",
      smoothness_orders = 1,
      max_degree = 2,
      num_knots = REFERENCE_DEFAULT_NUM_KNOTS,
      cv_folds = as.integer(config$cv_folds),
      cv_nlambda = as.integer(config$cv_nlambda),
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-4,
      prediction_bounds = "default"
    ),
    package_versions = list(
      hal9001 = package_version_string("hal9001"),
      glmnet = package_version_string("glmnet"),
      Matrix = package_version_string("Matrix"),
      dplyr = package_version_string("dplyr"),
      tidyr = package_version_string("tidyr")
    ),
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform
  )

  hash_r_object(manifest)
}


#' Return the checkpoint path for one run and signature.
#' @noRd
make_run_checkpoint_file <- function(checkpoint_dir, run_index,
                                     run_signature, attempt_id = NULL) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- if (is.null(attempt_id) || !nzchar(attempt_id)) {
    ""
  } else {
    paste0("_", safe_path_component(attempt_id))
  }
  file.path(
    checkpoint_dir,
    sprintf(
      "run_%04d_%s%s.rds",
      as.integer(run_index),
      substr(run_signature, 1L, 16L),
      suffix
    )
  )
}


#' Create a unique, readable attempt identifier.
make_run_attempt_id <- function() {
  safe_path_component(paste0(
    format(Sys.time(), "%Y%m%dT%H%M%OS6"),
    "-pid", Sys.getpid()
  ))
}


#' Read a valid, successful run checkpoint.
#' @noRd
read_successful_run_checkpoint <- function(path, expected_signature) {
  recover_atomic_target(path)
  if (!file.exists(path)) {
    return(NULL)
  }

  checkpoint <- tryCatch(readRDS(path), error = function(e) e)
  structurally_valid <- !inherits(checkpoint, "error") &&
    is.list(checkpoint) &&
    identical(checkpoint$schema, RUN_CHECKPOINT_SCHEMA) &&
    identical(checkpoint$run_signature, expected_signature) &&
    is.data.frame(checkpoint$timing) &&
    is.data.frame(checkpoint$fit_timings) &&
    is.data.frame(checkpoint$results)

  if (!isTRUE(structurally_valid)) {
    quarantine_invalid_file(path, "invalid-run-checkpoint")
    return(NULL)
  }

  complete <- nrow(checkpoint$timing) == 1L &&
    identical(as.character(checkpoint$status), "ok") &&
    identical(as.character(checkpoint$timing$status[[1]]), "ok") &&
    nrow(checkpoint$fit_timings) == 4L &&
    "fit_label" %in% names(checkpoint$fit_timings) &&
    setequal(
      as.character(checkpoint$fit_timings$fit_label),
      EXPECTED_HAL_FIT_LABELS
    ) &&
    nrow(checkpoint$results) == 2L &&
    all(c("Estimator", "Version") %in% names(checkpoint$results)) &&
    identical(
      sort(as.character(checkpoint$results$Version)),
      c("Naive", "Proposed")
    ) &&
    all(as.character(checkpoint$results$Estimator) == "AIPW")
  if (!isTRUE(complete)) {
    return(NULL)
  }

  checkpoint
}


#' Find the newest successful immutable checkpoint for one run.
#' @noRd
find_successful_run_checkpoint <- function(checkpoint_dir, run_index,
                                           run_signature) {
  canonical <- make_run_checkpoint_file(
    checkpoint_dir,
    run_index,
    run_signature
  )
  stem <- sub("[.]rds$", "", basename(canonical))
  candidates <- Sys.glob(file.path(checkpoint_dir, paste0(stem, "*.rds")))
  if (length(candidates) == 0L) {
    return(list(checkpoint = NULL, path = NULL))
  }

  info <- file.info(candidates)
  candidates <- candidates[order(info$mtime, decreasing = TRUE)]
  for (path in candidates) {
    checkpoint <- read_successful_run_checkpoint(path, run_signature)
    if (!is.null(checkpoint)) {
      return(list(checkpoint = checkpoint, path = path))
    }
  }
  list(checkpoint = NULL, path = NULL)
}


#' Atomically write per-run CSV shards from one checkpoint object.
#' @noRd
write_run_output_shards <- function(checkpoint, checkpoint_dir, run_index,
                                    run_signature, attempt_id = NULL) {
  shard_dir <- file.path(checkpoint_dir, "csv")
  dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(attempt_id) || !nzchar(attempt_id)) {
    attempt_id <- if (!is.null(checkpoint$attempt_id)) {
      checkpoint$attempt_id
    } else {
      "restored"
    }
  }
  stem <- sprintf(
    "run_%04d_%s_%s",
    as.integer(run_index),
    substr(run_signature, 1L, 16L),
    safe_path_component(attempt_id)
  )

  paths <- list(
    timing = file.path(shard_dir, paste0(stem, "_timing.csv")),
    fit_timing = file.path(shard_dir, paste0(stem, "_fit_timing.csv")),
    results = file.path(shard_dir, paste0(stem, "_results.csv"))
  )

  atomic_write_csv(checkpoint$timing, paths$timing, overwrite = FALSE)
  if (nrow(checkpoint$fit_timings) > 0L) {
    atomic_write_csv(
      checkpoint$fit_timings,
      paths$fit_timing,
      overwrite = FALSE
    )
  }
  if (nrow(checkpoint$results) > 0L) {
    atomic_write_csv(checkpoint$results, paths$results, overwrite = FALSE)
  }
  paths
}
