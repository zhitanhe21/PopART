test_that("Simulation 2 uses the global scheduler and creates figures", {
  project_dir <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  skip_if_not(file.exists(file.path(
    project_dir, "simulation", "sim_functions", "sim2.R"
  )))
  old_dir <- setwd(project_dir)
  on.exit(setwd(old_dir), add = TRUE)

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  source("simulation/sim_functions/global_scheduler.R")
  source("simulation/sim_functions/sim2.R")
  source("simulation/sim_analysis/sim2_analysis.R")

  study <- run_sim2(
    sample_sizes = data.frame(n_trial = 200L, n_auxiliary = 200L),
    run_id = "smoke_test",
    out_dir = file.path(tempdir(), "popart_simulation2"),
    mc_reps = 1L,
    base_seed = 100L,
    global_fit_slots = 2L,
    active_replicates = 1L,
    cv_lambdas = 10L
  )
  result <- study$results
  expect_equal(nrow(result), 2L)
  expect_equal(unique(result$Estimator), "AIPW")
  expect_setequal(result$Version, c("Naive", "Proposed"))

  sizes <- expand.grid(n_trial = c(500L, 5000L), n_aux = c(500L, 5000L))
  plot_results <- do.call(rbind, lapply(seq_len(nrow(sizes)), function(j) {
    do.call(rbind, lapply(seq_len(4L), function(i) {
      value <- result
      value$seed <- 10L * j + i
      value$n_trial <- sizes$n_trial[[j]]
      value$n_aux <- sizes$n_aux[[j]]
      value$etahat_0 <- value$etahat_0 + i / 1000
      value$etahat_1 <- value$etahat_1 - i / 1000
      value$rdhat <- value$etahat_1 - value$etahat_0
      value$rrhat <- value$etahat_1 / value$etahat_0
      value
    }))
  }))
  figures <- plot_sim2(
    plot_results, "smoke_test",
    file.path(tempdir(), "popart_simulation2_figures"), dpi = 72
  )
  expect_true(all(file.exists(figures)))
})
