#' aipwLocalDemo: Local AIPW and HAL Simulation Workflows
#'
#' Provides a reproducible local workflow for Simulation 2 augmented inverse
#' probability weighting (AIPW) estimators with highly adaptive lasso (HAL)
#' nuisance fits. The package includes laptop-aware resource planning, a
#' cross-run global HAL scheduler, recoverable checkpoints, and diagnostic
#' reporting helpers.
#'
#' The statistical data-generating process and estimator structure are derived
#' from upstream research code whose redistribution terms are still being
#' confirmed. See the package `LICENSE` file before sharing this package.
#'
#' @importFrom dplyr %>% any_of arrange bind_rows desc filter group_by matches
#'   mutate n rename select summarise
#' @importFrom ggplot2 aes annotate element_blank element_text facet_grid
#'   facet_wrap geom_abline geom_boxplot geom_hline geom_point ggplot ggsave
#'   labs position_jitter scale_color_manual scale_fill_manual
#'   scale_x_continuous scale_y_continuous theme theme_bw
#' @importFrom stats cov plogis predict qnorm rbinom rnorm var
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom utils flush.console read.csv
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "A", "C", "Estimator", "Param", "R", "Var", "Version", "W1", "W2",
  "W3", "X1", "Y", "Y0", "Y1", "cache_hit", "ci_cov", "ci_lower",
  "ci_upper", "clust", "cov_00", "cov_01", "cov_11", "emp_var", "est",
  "est_var", "eta_0", "eta_1", "etahat_0", "etahat_1", "fit_label",
  "max_effective_elapsed_sec", "mu0", "mu1", "n_aux", "n_both", "n_obs",
  "n_trial", "name", "pC", "pR", "param", "rd", "rdhat",
  "report_elapsed_sec", "rr", "rrhat", "run_index", "seed", "truth",
  "type", "value", "var_rd", "var_rr"
))
