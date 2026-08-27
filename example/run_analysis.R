# Import data ---------------------------------------------------------------

library(data.table)
library(popart)

example_dir <- dirname(normalizePath(sys.frame(1)$ofile))
data <- fread(file.path(example_dir, "popart_example.csv"), na.strings = "")


# Split trial and auxiliary samples ----------------------------------------

trial <- as.data.frame(data[S == 1])
auxiliary <- as.data.frame(data[S == 0])


# Fit PopART ---------------------------------------------------------------

fit <- popart::fit_popart(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "Y",
  treatment = "A",
  response = "R",
  censoring = "C",
  cluster = "cluster",
  covariates = c("X1", "X2", "W1", "W2")
)


# Results ------------------------------------------------------------------

print(fit)
print(summary(fit))
