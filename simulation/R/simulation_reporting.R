###############################################################################
###############################################################################

# Reporting helpers for the PopART published-reference study

###############################################################################
###############################################################################


#' Write a Markdown report for a Monte Carlo study.
#'
#' Reads the timing, fit timing, and estimator result CSVs written by
#' `run_monte_carlo_study()` and writes a compact Markdown report containing runtime,
#' fit timing, and estimate summaries.
#'
#' @param run_id Character run identifier used in output file names.
#' @param out_dir Character directory containing result and timing CSV files.
#'
#' @return Invisibly returns the path to the written Markdown report.
write_monte_carlo_report <- function(
    run_id,
    out_dir = file.path(getwd(), "simulation", "results")) {

  load_reference_study_packages()

  if (missing(run_id) || is.null(run_id) || !nzchar(run_id)) {
    stop("run_id is required", call. = FALSE)
  }

  prefix <- paste("monte_carlo", run_id, sep = "_")
  timing_file <- file.path(out_dir, paste0(prefix, "_timing.csv"))
  fit_timing_file <- file.path(out_dir, paste0(prefix, "_fit_timing.csv"))
  result_file <- file.path(out_dir, paste0(prefix, "_results.csv"))
  report_file <- file.path(out_dir, paste0(prefix, "_analysis_report.md"))

  timings <- if (file.exists(timing_file)) read.csv(timing_file, check.names = FALSE) else data.frame()
  fit_timings <- if (file.exists(fit_timing_file)) read.csv(fit_timing_file, check.names = FALSE) else data.frame()
  results <- if (file.exists(result_file)) read.csv(result_file, check.names = FALSE) else data.frame()

  runtime_display <- timings %>%
    select(any_of(c(
      "run_index", "n_trial", "n_aux", "status", "elapsed_min",
      "reserve_cores", "fit_workers", "cv_workers_per_fit", "backend",
      "use_cache", "resumed"
    )))

  fit_display <- data.frame()
  if (nrow(fit_timings) > 0L) {
    if ("effective_elapsed_sec" %in% names(fit_timings)) {
      fit_timings$report_elapsed_sec <- fit_timings$effective_elapsed_sec
    } else {
      fit_timings$report_elapsed_sec <- fit_timings$elapsed_sec
    }
    fit_display <- fit_timings %>%
      group_by(fit_label) %>%
      summarise(
        n_runs = n(),
        cache_hits = sum(cache_hit %in% TRUE, na.rm = TRUE),
        mean_effective_elapsed_sec = mean(
          report_elapsed_sec,
          na.rm = TRUE
        ),
        max_effective_elapsed_sec = max(
          report_elapsed_sec,
          na.rm = TRUE
        ),
        mean_n_obs = mean(n_obs, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(max_effective_elapsed_sec))
  }

  estimate_display <- data.frame()
  if (nrow(results) > 0L) {
    estimate_display <- results %>%
      mutate(
        rd_error = rdhat - rd,
        rr_error = rrhat - rr
      ) %>%
      select(any_of(c(
        "run_index", "n_trial", "n_aux", "Version",
        "rdhat", "rd", "rd_error", "rrhat", "rr", "rr_error"
      ))) %>%
      arrange(run_index, Version)
  }

  report <- c(
    "# PopART Published-Reference Monte Carlo Report",
    "",
    paste0("- Run ID: `", run_id, "`"),
    paste0("- Generated: `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "`"),
    paste0("- Output directory: `", out_dir, "`"),
    "",
    "## Runtime",
    "",
    md_table(runtime_display),
    "",
    "## Fit Timing",
    "",
    md_table(fit_display),
    "",
    "## Estimates",
    "",
    md_table(estimate_display, digits = 6),
    ""
  )

  writeLines(report, report_file)
  cat("Wrote report: ", report_file, "\n", sep = "")

  invisible(report_file)
}


#' Create a placeholder ggplot for unavailable diagnostics.
#'
#' @param title Character plot title.
#' @param message Character message shown in the plot body.
#'
#' @return A `ggplot` object.
make_placeholder_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = message,
      size = 4
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}


#' Summarize Monte Carlo estimates for plotting.
#'
#' Converts wide estimator output into long-form truth, estimate, and variance
#' summaries for eta(0), eta(1), RD, and RR.
#'
#' @param results Data frame returned by `run_monte_carlo_study()` or read from the
#'   corresponding results CSV.
#'
#' @return Data frame with bias, empirical variance, average estimated variance,
#'   MSE, and empirical CI coverage by estimator, sample size, and parameter.
summarize_monte_carlo <- function(results) {
  load_reference_study_packages()
  results %>%
    select(seed, Version, Estimator, n_trial, n_aux,
           eta_0, eta_1, rd, rr,
           etahat_0, etahat_1, rdhat, rrhat,
           cov_00, cov_01, cov_11, var_rd, var_rr) %>%
    rename(
      truth_eta0 = eta_0,
      truth_eta1 = eta_1,
      truth_rd = rd,
      truth_rr = rr,
      est_eta0 = etahat_0,
      est_eta1 = etahat_1,
      est_rd = rdhat,
      est_rr = rrhat,
      Var_eta0 = cov_00,
      Var_eta1 = cov_11,
      Var_rd = var_rd,
      Var_rr = var_rr
    ) %>%
    select(-cov_01) %>%
    tidyr::pivot_longer(
      cols = matches("^(truth|est|Var)_"),
      names_to = c("type", "param"),
      names_sep = "_",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(
      names_from = type,
      values_from = value,
      id_cols = c(seed, Version, Estimator, n_trial, n_aux, param)
    ) %>%
    mutate(
      ci_lower = est - qnorm(0.975) * sqrt(Var),
      ci_upper = est + qnorm(0.975) * sqrt(Var)
    ) %>%
    group_by(Version, Estimator, n_trial, n_aux, param) %>%
    summarise(
      bias = mean(est - truth, na.rm = TRUE),
      emp_var = var(est, na.rm = TRUE),
      est_var = mean(Var, na.rm = TRUE),
      mse = mean((est - truth)^2, na.rm = TRUE),
      ci_cov = mean(truth >= ci_lower & truth <= ci_upper, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Param = factor(
        param,
        levels = c("eta0", "eta1", "rd", "rr"),
        labels = c("eta(0)", "eta(1)", "RD", "RR")
      )
    )
}


#' Write the three standard reference-study figures.
#'
#' Generates RD estimates, variance calibration, and confidence-interval
#' coverage figures. When too few seeds are available to compute empirical
#' variance, a placeholder variance figure is written so the workflow still
#' produces the expected three files.
#'
#' @param run_id Character run identifier used in output file names.
#' @param out_dir Character directory containing the results CSV.
#' @param figure_dir Character directory where figures are written.
#' @param width,height Numeric figure dimensions in inches.
#' @param dpi Integer output image resolution.
#'
#' @return Invisibly returns a named character vector of written figure paths.
write_monte_carlo_figures <- function(
    run_id,
    out_dir = file.path(getwd(), "simulation", "results"),
    figure_dir = file.path(out_dir, "figures"),
    width = 8,
    height = 6,
    dpi = 600) {

  load_reference_study_packages()

  if (missing(run_id) || is.null(run_id) || !nzchar(run_id)) {
    stop("run_id is required", call. = FALSE)
  }

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- paste("monte_carlo", run_id, sep = "_")
  result_file <- file.path(out_dir, paste0(prefix, "_results.csv"))
  if (!file.exists(result_file)) {
    stop("Could not find result file: ", result_file, call. = FALSE)
  }

  sim_res <- read.csv(result_file, check.names = FALSE) %>%
    mutate(
      Estimator = factor(Estimator, levels = c("AIPW")),
      Version = factor(Version, levels = c("Naive", "Proposed")),
      Trial_Size = factor(n_trial),
      Aux_Size = factor(n_aux),
      Size = factor(paste0(n_trial, "_", n_aux))
    )

  pal <- c("#FF6800", "#803E75", "#C10020", "#FFB300")
  rd <- mean(sim_res$rd, na.rm = TRUE)

  estimate_plot <- sim_res %>%
    ggplot(aes(
      x = interaction(n_aux, n_trial),
      y = rdhat,
      color = Version,
      fill = Version
    )) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(
      position = position_jitter(width = 0.08, height = 0),
      alpha = 0.65,
      size = 1.6
    ) +
    geom_hline(yintercept = rd, linetype = "dashed") +
    facet_wrap(~ Version) +
    scale_color_manual(values = pal[seq_len(nlevels(sim_res$Version))]) +
    scale_fill_manual(values = pal[seq_len(nlevels(sim_res$Version))]) +
    labs(
      y = expression(hat(RD)),
      x = "Auxiliary and Trial Sample Sizes"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

  summary_table <- summarize_monte_carlo(sim_res)

  variance_data <- summary_table %>%
    filter(is.finite(emp_var), is.finite(est_var), emp_var > 0, est_var > 0) %>%
    mutate(n_both = factor(paste0(n_aux, "_", n_trial)))

  if (nrow(variance_data) > 0L) {
    variance_plot <- variance_data %>%
      ggplot(aes(
        x = emp_var,
        y = est_var,
        color = n_both,
        shape = Param
      )) +
      scale_x_continuous(transform = "log10") +
      scale_y_continuous(transform = "log10") +
      scale_color_manual(values = pal[seq_len(nlevels(variance_data$n_both))]) +
      geom_abline(linetype = "dashed") +
      geom_point(size = 3) +
      facet_grid(Estimator ~ Version) +
      labs(
        y = "Average Estimated Variance",
        x = "Empirical Variance",
        color = "Auxiliary and Trial Sample Size",
        shape = "Parameter"
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        legend.box = "vertical"
      )
  } else {
    variance_plot <- make_placeholder_plot(
      "Average Estimated Variance vs Empirical Variance",
      "Empirical variance needs at least two simulation seeds per group."
    )
  }

  coverage_data <- summary_table %>%
    filter(is.finite(ci_cov))

  if (nrow(coverage_data) > 0L) {
    coverage_plot <- coverage_data %>%
      ggplot(aes(
        x = interaction(n_aux, n_trial),
        y = ci_cov,
        color = Param,
        shape = Param
      )) +
      geom_point(size = 3) +
      geom_hline(yintercept = 0.95, linetype = "dashed") +
      scale_y_continuous(limits = c(0, 1)) +
      scale_color_manual(values = pal[seq_len(nlevels(coverage_data$Param))]) +
      facet_grid(Estimator ~ Version) +
      labs(
        y = "Empirical CI Coverage",
        x = "Auxiliary and Trial Sample Sizes",
        color = "Estimand",
        shape = "Estimand"
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )
  } else {
    coverage_plot <- make_placeholder_plot(
      "Empirical CI Coverage",
      "Coverage could not be computed from the available result rows."
    )
  }

  estimate_file <- file.path(figure_dir, paste0(prefix, "_estimates.png"))
  variance_file <- file.path(figure_dir, paste0(prefix, "_variance.png"))
  coverage_file <- file.path(figure_dir, paste0(prefix, "_confidence.png"))

  ggsave(estimate_file, estimate_plot, dpi = dpi, width = width, height = height)
  ggsave(variance_file, variance_plot, dpi = dpi, width = width, height = height)
  ggsave(coverage_file, coverage_plot, dpi = dpi, width = width, height = height)

  cat("Wrote estimate figure: ", estimate_file, "\n", sep = "")
  cat("Wrote variance figure: ", variance_file, "\n", sep = "")
  cat("Wrote confidence figure: ", coverage_file, "\n", sep = "")

  invisible(c(
    estimates = estimate_file,
    variance = variance_file,
    confidence = coverage_file
  ))
}
