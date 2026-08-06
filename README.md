# popart

`popart` analyzes a randomized-trial sample together with an auxiliary sample
representing the target population. The package fits highly adaptive lasso
(HAL) nuisance regressions and reports two augmented inverse probability
weighting (AIPW) analyses: a trial-only comparator and an analysis that also
uses the auxiliary sample.

The installed package analyzes data supplied by the user. It does not generate
research data or run a Monte Carlo experiment. Repository-only simulation and
method-validation tools are described in the
[simulation workflow](#repository-simulation-workflow) below.

## Installation

From the package directory:

```sh
R CMD INSTALL .
```

The source currently uses the repository's existing license notice. Confirm
the terms in [`LICENSE`](LICENSE) before redistributing this package or making
the repository public.

## Data layout

`fit_popart()` takes two data frames:

- `trial_data` contains the outcome, randomized treatment, response indicator,
  censoring indicator, and baseline covariates;
- `auxiliary_data` contains the same baseline covariates and may contain an
  auxiliary sampling-weight column. It does not need treatment, response,
  censoring, or outcome columns.

Treatment values are supplied in control-then-treated order. Response and
censoring indicators use 1 for "responded" and "censored", respectively. The
outcome, response, and censoring variables must be binary. Covariates must be
numeric or logical and must not be missing; categorical variables should be
encoded before fitting.

The package accepts data frames, not file paths. For example, CSV files can be
read with `read.csv()` before calling the estimator.

## Basic analysis

The package includes two small synthetic CSV files solely to demonstrate the
input structure:

```r
library(popart)

trial <- read.csv(
  system.file("extdata", "example_trial.csv", package = "popart")
)
auxiliary <- read.csv(
  system.file("extdata", "example_auxiliary.csv", package = "popart")
)

fit <- fit_popart(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "outcome",
  treatment = "treatment",
  response = "responded",
  censoring = "censored",
  covariates = c(
    "baseline_risk",
    "baseline_binary",
    "baseline_score_1",
    "baseline_score_2"
  ),
  auxiliary_weight = "survey_weight",
  treatment_values = c(0, 1)
)

summary(fit)
coef(fit, estimator = "trial_auxiliary")
vcov(fit, estimator = "trial_auxiliary")
confint(fit, estimator = "trial_auxiliary")
```

The returned `popart_fit` object contains estimates, covariance matrices,
nuisance-fit timing, prediction diagnostics, and the analysis specification.
Fitted nuisance models are omitted by default to keep the result smaller and
avoid retaining input-derived design matrices. Set `keep_nuisance_fits = TRUE`
in `popart_control()` only when those models are needed after estimation.

The included CSV files are a fixed synthetic snapshot generated from the
repository's reference data-generating process with seed `20260805`, 20
clusters, 200 trial records, and 200 auxiliary records. Missing trial outcomes
and censoring indicators are represented by blank CSV fields when they are not
observed. The auxiliary file intentionally omits trial-only variables.

## Computation controls

The default is serial and is suitable as a conservative starting point. HAL
fits can instead be parallelized at one of two levels:

```r
# Run two complete nuisance fits concurrently.
control <- popart_control(n_fit_workers = 2, n_cv_workers = 1)

# Or run complete fits serially and parallelize their CV folds.
control <- popart_control(n_fit_workers = 1, n_cv_workers = 3)
```

Do not enable both levels at once. Worker counts refer to operating-system
workers, while the relationship between those workers and physical CPU cores
depends on the computer and operating system. On a laptop, begin with one or
two workers and leave capacity for the operating system and other programs.

The tuning defaults are three CV folds and `num_knots = c(5, 3)`. These are
computational defaults, not universally optimal scientific choices. Record any
changes to tuning parameters as part of the analysis specification.

## Repository simulation workflow

The repository-level [`simulation/`](simulation/) directory is excluded from
the installed package. It contains research tools used to generate reference
data, run Monte Carlo studies, validate the global HAL scheduler, and recreate
simulation reports. Applied analyses do not need this directory.

Its main components are:

- `simulation/R/`: data generation, reference estimators, global scheduling,
  persistence, and reporting functions;
- `simulation/scripts/run_monte_carlo.R`: command-line Monte Carlo entry point;
- `simulation/scripts/analyze_results.R`: regenerate reports and figures;
- `simulation/scripts/create_example_data.R`: recreate the fixed package CSVs;
- `simulation/tests/formal_api_equivalence.R`: compare `fit_popart()` with the
  preserved reference equations;
- `simulation/results/`: ignored generated results, checkpoints, and caches.

### Generate, export, and analyze reference data

The following repository-only example shows how simulated data pass through
the same two-data-frame API used for real imported data:

```r
library(popart)

source("simulation/R/load_simulation.R")
study <- new.env(parent = globalenv())
source_reference_study(".", envir = study)

generated <- study$generate_reference_data(
  m = 20,
  n_trial = 500,
  n_auxiliary = 500,
  p_resp = 0.5,
  p_cens = 0.3,
  seed = 20260805
)

covariates <- c("X1", "W1", "W2", "W3")
trial_generated <- generated$dat[generated$dat$S == 1, ]
auxiliary_generated <- generated$dat[
  generated$dat$S == 0,
  c(covariates, "wt")
]

teaching_dir <- file.path(tempdir(), "popart-teaching-data")
dir.create(teaching_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  trial_generated,
  file.path(teaching_dir, "trial.csv"),
  row.names = FALSE
)
write.csv(
  auxiliary_generated,
  file.path(teaching_dir, "auxiliary.csv"),
  row.names = FALSE
)

trial <- read.csv(file.path(teaching_dir, "trial.csv"))
auxiliary <- read.csv(file.path(teaching_dir, "auxiliary.csv"))

fit <- fit_popart(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "Y",
  treatment = "A",
  response = "R",
  censoring = "C",
  covariates = covariates,
  auxiliary_weight = "wt",
  control = popart_control(random_seed = 20260805)
)

summary(fit)
confint(fit)
```

The generator is used only to create teaching or validation data. With a real
study, start at the `read.csv()` step and replace the paths and column names.

### Run the Monte Carlo study

Sample sizes are supplied explicitly; there are no `small`, `medium`, or
`large` workload labels:

```text
Rscript simulation/scripts/run_monte_carlo.R --n-trial 500,5000 --n-auxiliary 500,5000 --replicates 20 --run-id reference_20
```

Inspect the scheduler plan without starting HAL fits:

```text
Rscript simulation/scripts/run_monte_carlo.R --n-trial 500,5000 --n-auxiliary 500,5000 --replicates 20 --dry-run
```

After installing `popart`, verify that the formal API still reproduces the
preserved reference implementation:

```sh
Rscript simulation/tests/formal_api_equivalence.R
```

The synthetic CSV files in `inst/extdata/` are teaching fixtures, not benchmark
results or real participant data. Package users can also see
[`vignettes/getting-started.Rmd`](vignettes/getting-started.Rmd) for the shorter
input-to-estimate walkthrough installed with the package.
