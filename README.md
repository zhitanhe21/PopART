# popart: Causal Inference for Cluster-Randomized Trials with Differential Nonresponse

## Background

`popart` estimates population-level causal effects from a cluster-randomized
trial with differential individual nonresponse and outcome censoring.

| Design feature | Statistical role |
|---|---|
| Cluster-level treatment | Every individual in a randomized cluster receives the same treatment assignment |
| Individual response | Only some invited individuals respond and provide individual-level data |
| Outcome censoring | A responder's outcome may remain unobserved |
| Auxiliary sample | Supplies the baseline covariate distribution of the target population |

If response differs by treatment or baseline covariates, complete trial cases
may not represent the population in the randomized clusters. The proposed
estimators combine observed trial outcomes with the auxiliary covariates to
address this selection problem.

For binary outcomes, the causal targets are:

| Parameter | Definition | Meaning |
|---|---|---|
| `eta(0)` | `E[Y(0)]` | Population risk if all clusters received control |
| `eta(1)` | `E[Y(1)]` | Population risk if all clusters received treatment |
| `RD` | `eta(1) - eta(0)` | Causal risk difference |
| `RR` | `eta(1) / eta(0)` | Causal risk ratio |

## Data structure

The package receives trial and auxiliary data frames separately. The example
CSV uses `S` to identify the two samples before splitting them.

| Variable | Level | Meaning | Trial | Auxiliary |
|---|---|---|---:|---:|
| `id` | Individual | Observation identifier | yes | yes |
| `cluster` | Cluster | Randomization unit | yes | yes |
| `S` | Individual | `1` = trial, `0` = auxiliary; used only to split the CSV | yes | yes |
| `A` | Cluster | Treatment assignment; constant within cluster | yes | not required |
| `R` | Individual | `1` = trial respondent | yes | no |
| `C` | Individual | `1` = censored outcome among respondents | yes | no |
| `Y` | Individual | Binary outcome, observed when `R = 1` and `C = 0` | yes | no |
| `X1`, `X2` | Cluster | Cluster-level baseline covariates | yes | yes |
| `W1`, `W2` | Individual | Individual-level baseline covariates | column required; responders used | yes |

**Required cluster structure:** within `trial_data`, every row with the same
`cluster` value must have the same non-missing treatment `A`. The auxiliary
sample does not need a treatment column. Data with treatment varying between
individuals in the same cluster do not match this cluster-randomized design.

Causal interpretation requires cluster randomization, positivity, sufficient
measured covariates for response and censoring, and a representative auxiliary
sample or valid auxiliary sampling weights.

### Auxiliary sampling weights

If the auxiliary survey supplies sampling weights, pass that column through
`auxiliary_weight`. Otherwise leave it as `NULL`; the package does not estimate
or invent sampling weights.

## Installation

Run these commands from the `PopART_full` folder:

```r
install.packages(c("hal9001", "numDeriv", "data.table"))
install.packages(".", repos = NULL, type = "source")
library(popart)
```

The simulation scripts additionally use:

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggh4x", "legendry", "pkgload"
))
```

## Quick start

Open `example/run_analysis.R` and click **Run**, or run:

```r
source("example/run_analysis.R")
```

The complete analysis is:

```r
library(data.table)
library(popart)

data <- fread("example/popart_example.csv", na.strings = "")
trial <- as.data.frame(data[S == 1])
auxiliary <- as.data.frame(data[S == 0])

fit <- fit_popart(
  trial_data = trial,
  auxiliary_data = auxiliary,
  outcome = "Y",
  treatment = "A",
  response = "R",
  censoring = "C",
  cluster = "cluster",
  covariates = c("X1", "X2", "W1", "W2")
)

print(fit)
summary(fit)
```

### Example data

The example contains 1,000,000 rows from 200 clusters: 500,000 trial rows and
500,000 auxiliary rows. It is analysis-ready and is split by `S` before
`fit_popart()` is called. `NA` marks unavailable response, censoring, and
outcome fields.

![Trial and auxiliary rows from the example data](example/figures/example_data_structure.png)

### Example output

`print(fit)` gives a concise comparison, while `summary(fit)` reports the
parameter-level estimates, standard errors, and 95% confidence intervals.
The upper table is the quick comparison and the lower table is the detailed
inference output.

![PopART example estimates and confidence intervals](example/figures/example_output.png)

## Main function and parameters

Only one function is needed for a real-data analysis.

| Function | Purpose | Output |
|---|---|---|
| `fit_popart()` | Fit Naive and Proposed G-formula, IPW, and AIPW estimators | `estimates`, `variance`, and `covariance` |
| `print(fit)` | Display the six estimator/version combinations | `eta(0)`, `eta(1)`, `RD`, and `RR` |
| `summary(fit)` | Display all parameter-level results | Estimate, standard error, and 95% CI |

### Parameters of `fit_popart()`

| Parameter | Meaning | Default |
|---|---|---|
| `trial_data` | Trial observations | required |
| `auxiliary_data` | Auxiliary target-population observations | required |
| `outcome` | Binary outcome column | required |
| `treatment` | Cluster-level treatment column | required |
| `response` | Trial response indicator column | required |
| `censoring` | Outcome-censoring indicator column | required |
| `cluster` | Cluster identifier column | required |
| `covariates` | Shared cluster- and individual-level baseline covariates | required |
| `auxiliary_weight` | Optional survey-weight column supplied with the auxiliary data | `NULL` |
| `treatment_values` | Control and treatment values | `c(0, 1)` |
| `n_cv_folds` | Number of HAL cross-validation folds | `5` |
| `n_lambda_values` | Number of HAL lasso penalties | `50` |
| `random_seed` | Seed used to create cluster-grouped HAL folds | `1` |

## Simulation 1: parametric model specification

### Purpose

Simulation 1 compares Naive and Proposed G-formula, IPW, and AIPW estimators
when the parametric outcome and selection models are correctly specified or
misspecified.

### DGP and settings

| Component | Setting |
|---|---|
| Clusters | 20; half assigned to control and half to treatment |
| Trial sizes | 500 and 5000 |
| Auxiliary sizes | 500 and 5000 |
| Monte Carlo replicates | 20 per sample-size combination |
| Cluster covariate | `X1`: fixed cluster-level risk |
| Individual covariates | `W1 ~ Bernoulli(0.5)` in trial; `W1 ~ Bernoulli(0.75)` in auxiliary; `W2 ~ Normal(0,1)` |
| Response | Binary; depends on the treatment-by-`W1` interaction; marginal rate 0.5 |
| Censoring | Binary; depends on treatment and `W1`; rate among responders 0.3 |
| Control outcome | `logit P(Y(0)=1) = -1 + 2W1 + 0.5W2 + 0.25X1` |
| Treated outcome | `logit P(Y(1)=1) = -W1 - 0.5W2` |

| Scenario | Correct model | Misspecified model |
|---|---|---|
| Outcome `mu` | Includes `X1`, `W1`, `W2`, treatment, and treatment interactions | Omits `W1` and its treatment interaction |
| Selection/censoring `pi` | Includes `X1`, `W1`, and `W2` | Omits `W1` |

Each replicate combines trial rows (`S = 1`) with auxiliary covariate rows
(`S = 0`). Auxiliary `R`, `C`, and `Y` values are internal placeholders and
are not treated as observed outcomes.

### Current results

**RD estimates.** The boxplots show the Monte Carlo distribution of estimated
risk differences. The dashed line is the true RD; gray panels mark scenarios
in which the corresponding estimator is not expected to be consistent.

![Simulation 1 estimates](simulation/sim_figures/sim1/simulation1_parametric_mc20_estimates.png)

**Variance calibration.** Each point compares empirical variance with average
estimated variance. Points close to the dashed diagonal indicate agreement.

![Simulation 1 variance](simulation/sim_figures/sim1/simulation1_parametric_mc20_variance.png)

**Confidence-interval coverage.** Points show empirical coverage for
`eta(0)`, `eta(1)`, RD, and RR; the dashed line marks the nominal 95% level.

![Simulation 1 confidence interval coverage](simulation/sim_figures/sim1/simulation1_parametric_mc20_confidence.png)

### Reproduction

Open `simulation/sim_scripts/ss1.R` and click **Run**. The script generates the
data, fits the estimators, and writes the result table to
`simulation/sim_data/sim1/` and the figures to
`simulation/sim_figures/sim1/`.

## Simulation 2: HAL AIPW with the global scheduler

### Purpose

Simulation 2 compares Naive and Proposed AIPW estimators when nonlinear
response, censoring, and outcome functions are learned with HAL. It also tests
the global scheduler used to distribute HAL fits across Monte Carlo replicates.

### DGP and settings

| Component | Setting |
|---|---|
| Clusters | 20; half assigned to control and half to treatment |
| Trial sizes | 500 and 5000 |
| Auxiliary sizes | 500 and 5000 |
| Monte Carlo replicates | 20 per sample-size combination |
| Cluster covariate | `X1`: fixed cluster-level risk |
| Individual covariates | `W1` binary; `W2` and `W3` standard normal |
| Response | Nonlinear in treatment, `W1`, `W2`, and `W3`; marginal rate 0.5 |
| Censoring | Nonlinear in treatment, `W1`, `W2`, and `W3`; rate among responders 0.3 |
| Control outcome | Nonlinear risk depending on `W1` |
| Treated outcome | Nonlinear risk depending on `W1`, `sin(W2)`, `W3^2`, and `X1` |
| Estimators | Naive HAL AIPW and Proposed HAL AIPW |

Simulation 2 uses the same trial/auxiliary layout as Simulation 1, adds `W3`,
and fits HAL with `X1`, `W1`, `W2`, and `W3`.

### Global scheduler

Simulation 2 automatically uses `detectCores(logical = TRUE) - 2` workers in
one global pool and processes Monte Carlo replicates in batches.

### Current results

**RD estimates.** The boxplots compare Naive and Proposed HAL AIPW across the
four sample-size combinations. The dashed line is the true RD.

![Simulation 2 estimates](simulation/sim_figures/sim2/monte_carlo_hal_aipw_mc20_estimates.png)

**Variance calibration.** Each point compares empirical variance with average
estimated variance for the four causal parameters. The dashed diagonal marks
perfect agreement.

![Simulation 2 variance](simulation/sim_figures/sim2/monte_carlo_hal_aipw_mc20_variance.png)

**Confidence-interval coverage.** Points show empirical 95% CI coverage by
sample size and parameter; the dashed line marks the nominal 0.95 level.

![Simulation 2 confidence interval coverage](simulation/sim_figures/sim2/monte_carlo_hal_aipw_mc20_confidence.png)

### Reproduction

Open `simulation/sim_scripts/ss2.R` and click **Run**. CPU allocation, HAL
fitting, and output generation are handled by the script. Results are written
to `simulation/sim_data/sim2/` and figures to `simulation/sim_figures/sim2/`.
