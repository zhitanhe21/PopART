#!/usr/bin/env Rscript

###############################################################################
# Simulation 2: global-scheduler HAL AIPW Naive/Proposed comparison
# Click Run to execute the complete simulation.
###############################################################################

project_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", ".."))
setwd(project_dir)

# Packages -------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)


# Code -----------------------------------------------------------------------

pkgload::load_all(project_dir)
source("simulation/sim_functions/global_scheduler.R")
source("simulation/sim_functions/sim2.R")
source("simulation/sim_analysis/sim2_analysis.R")


# Parameters -----------------------------------------------------------------

sample_sizes <- expand.grid(
  n_trial = c(500L, 5000L),
  n_auxiliary = c(500L, 5000L)
)
mc_reps <- 20L
global_fit_slots <- parallel::detectCores(logical = TRUE) - 2L
active_replicates <- ceiling(global_fit_slots / 2L)

run_id <- "hal_aipw_mc20"
data_dir <- file.path(project_dir, "simulation", "sim_data", "sim2")
figure_dir <- file.path(project_dir, "simulation", "sim_figures", "sim2")


# Pipeline -------------------------------------------------------------------

simulation <- run_sim2(
  sample_sizes = sample_sizes,
  run_id = run_id,
  out_dir = data_dir,
  mc_reps = mc_reps,
  global_fit_slots = global_fit_slots,
  active_replicates = active_replicates
)

figures <- plot_sim2(
  simulation$results,
  run_id,
  figure_dir
)

cat(
  "\nCompleted", mc_reps, "Monte Carlo replicates for each of",
  nrow(sample_sizes), "sample-size combinations\n"
)
cat("Results:", simulation$result_file, "\n")
cat("Figures:\n")
print(figures)
