# Simulation 1 estimates, variances, and confidence coverage -------------

# 3. Summarize and plot the simulation results.
plot_sim1 <- function(results, run_id, figure_dir, dpi = 600) {
  summary <- results %>%
    dplyr::mutate(
      Mu = dplyr::if_else(
        Version == "Naive", "-",
        dplyr::if_else(mu_correct, "Correct Mu", "Incorrect Mu")
      ),
      Pi = dplyr::if_else(
        Version == "Naive", "-",
        dplyr::if_else(pi_correct, "Correct Pi", "Incorrect Pi")
      )
    ) %>%
    dplyr::filter(Version == "Proposed" | (mu_correct & pi_correct)) %>%
    dplyr::select(
      seed, Version, Estimator, Mu, Pi, n_trial, n_aux,
      eta_0, eta_1, rd, rr, etahat_0, etahat_1, rdhat, rrhat,
      cov_00, cov_11, var_rd, var_rr
    ) %>%
    dplyr::rename(
      truth_eta0 = eta_0, truth_eta1 = eta_1,
      truth_rd = rd, truth_rr = rr,
      est_eta0 = etahat_0, est_eta1 = etahat_1,
      est_rd = rdhat, est_rr = rrhat,
      Var_eta0 = cov_00, Var_eta1 = cov_11,
      Var_rd = var_rd, Var_rr = var_rr
    ) %>%
    tidyr::pivot_longer(
      dplyr::matches("^(truth|est|Var)_"),
      names_to = c("type", "parameter"), names_sep = "_"
    ) %>%
    tidyr::pivot_wider(names_from = type, values_from = value) %>%
    dplyr::mutate(
      lower = est - 1.96 * sqrt(Var),
      upper = est + 1.96 * sqrt(Var)
    ) %>%
    dplyr::group_by(Version, Estimator, Mu, Pi, n_trial, n_aux, parameter) %>%
    dplyr::summarise(
      bias = mean(est - truth),
      empirical_variance = stats::var(est),
      estimated_variance = mean(Var),
      mse = mean((est - truth)^2),
      coverage = mean(truth >= lower & truth <= upper),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW")),
      Version = factor(Version, c("Naive", "Proposed")),
      Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
      Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu")),
      parameter = factor(
        parameter,
        c("eta0", "eta1", "rd", "rr"),
        c("eta(0)", "eta(1)", "RD", "RR")
      )
    )

  sample_size <- function(auxiliary, trial) {
    interaction(
      factor(auxiliary, levels = sort(unique(auxiliary))),
      factor(trial, levels = sort(unique(trial))),
      sep = ".", drop = TRUE
    )
  }
  formatted <- results %>%
    dplyr::mutate(
      Estimator = factor(Estimator, c("G-Formula", "IPW", "AIPW")),
      Version = factor(Version, c("Naive", "Proposed")),
      Mu = dplyr::if_else(
        Version == "Naive", "-",
        dplyr::if_else(mu_correct, "Correct Mu", "Incorrect Mu")
      ),
      Pi = dplyr::if_else(
        Version == "Naive", "-",
        dplyr::if_else(pi_correct, "Correct Pi", "Incorrect Pi")
      ),
      consistent = Version == "Proposed" & (
        (Estimator == "G-Formula" & mu_correct) |
        (Estimator == "IPW" & pi_correct) |
        (Estimator == "AIPW" & (mu_correct | pi_correct))
      ),
      Sample_Size = sample_size(n_aux, n_trial)
    ) %>%
    dplyr::filter(Version == "Proposed" | (mu_correct & pi_correct)) %>%
    dplyr::mutate(
      Pi = factor(Pi, c("-", "Incorrect Pi", "Correct Pi")),
      Mu = factor(Mu, c("-", "Incorrect Mu", "Correct Mu"))
    )
  colors <- c("G-Formula" = "#FF6800", "IPW" = "#803E75", "AIPW" = "#C10020")
  shade <- formatted %>%
    dplyr::distinct(Estimator, Version, Pi, Mu, consistent) %>%
    dplyr::filter(!consistent) %>%
    dplyr::mutate(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

  estimates <- ggplot2::ggplot(
    formatted,
    ggplot2::aes(Sample_Size, rdhat, color = Estimator, fill = Estimator)
  ) +
    ggplot2::geom_rect(
      data = shade,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#e5e4e2"
    ) +
    ggplot2::geom_boxplot(alpha = 0.5) +
    ggplot2::geom_hline(yintercept = mean(formatted$rd), linetype = "dashed") +
    ggh4x::facet_nested(Estimator ~ Version + Pi + Mu) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "Auxiliary and Trial Sample Sizes", y = expression(hat(RD))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none") +
    ggplot2::guides(x = legendry::guide_axis_nested())

  variance_data <- summary %>%
    dplyr::filter(
      n_trial == n_aux,
      empirical_variance > 0,
      estimated_variance > 0
    )
  variance <- ggplot2::ggplot(
    variance_data,
    ggplot2::aes(
      empirical_variance, estimated_variance,
      color = factor(n_trial), shape = parameter
    )
  ) +
    ggplot2::geom_abline(linetype = "dashed") +
    ggplot2::geom_point(size = 3) +
    ggh4x::facet_nested(Estimator ~ Version + Pi + Mu) +
    ggplot2::scale_x_continuous(
      transform = "log10",
      breaks = c(0.001, 0.01),
      labels = c("0.001", "0.01")
    ) +
    ggplot2::scale_y_continuous(transform = "log10") +
    ggplot2::labs(
      x = "Empirical Variance", y = "Average Estimated Variance",
      color = "Trial and Auxiliary Sample Size", shape = "Parameter"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")

  coverage_data <- summary %>%
    dplyr::mutate(Sample_Size = sample_size(n_aux, n_trial))
  coverage <- ggplot2::ggplot(
    coverage_data,
    ggplot2::aes(Sample_Size, coverage, color = parameter, shape = parameter)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed") +
    ggh4x::facet_nested(Estimator ~ Version + Pi + Mu) +
    ggplot2::labs(
      x = "Auxiliary and Trial Sample Sizes", y = "Empirical CI Coverage",
      color = "Parameter", shape = "Parameter"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom") +
    ggplot2::guides(x = legendry::guide_axis_nested())

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(
    figure_dir,
    paste0("simulation1_", run_id, c("_estimates.png", "_variance.png", "_confidence.png"))
  )
  ggplot2::ggsave(files[1], estimates, width = 8, height = 6, dpi = dpi)
  ggplot2::ggsave(files[2], variance, width = 10, height = 6, dpi = dpi)
  ggplot2::ggsave(files[3], coverage, width = 8, height = 6, dpi = dpi)
  invisible(files)
}
