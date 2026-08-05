# Provenance

This version-2 local demo is derived from `aipw_local_demo` and the optimized
Fit+CV reuse Simulation 2 code in the HAL workspace. It is intentionally stored
in `aipw_local_demo_2` so the original demo, benchmark, and Slurm-oriented
workflows remain unchanged.

## Preserved DGP and Estimator Structure

The demo preserves the following core pieces from the optimized workflow:

- Simulation 2 clustered data generation.
- Shared outcome HAL fit reused by both naive and proposed AIPW estimators.
- Proposed Q regressions for the augmentation terms.
- Naive censoring regression for the baseline estimator.
- HAL still uses `smoothness_orders = 1`, `max_degree = 2`, and glmnet. The
  current experimental laptop defaults are `num_knots = c(5, 3)` and
  `cv_folds = 3`; Brian's original settings were `c(25, 10)` and four folds.
  These reductions change the HAL dictionary and CV-selected fit, so current
  results are not expected to be numerically identical to Brian's results.
  The programmatic API retains `cv_nlambda = 50`; the direct CLI intentionally
  defaults to `cv_nlambda = 5` for a quick local check.
- Fixed fold IDs for reproducible HAL CV.
- Fixed Simulation 2 response/censoring targets `p_resp = 0.5` and
  `p_cens = 0.3`, matching the intercepts in the original data-generating
  mechanism.
- Result assembly for `eta_0`, `eta_1`, RD, RR, covariance terms, estimated
  variances, and truth values.

## Local Demo Adaptations

The demo adds local-run convenience around the same statistical core:

- Brian-style organization under `simulation/sim_scripts/` and
  `simulation/sim_analysis/`,
  `simulation/sim_data/`, and `simulation/sim_figures/`.
- Package-loadable workflow and report functions under `R/`.
- A direct script, `simulation/sim_scripts/ss2_demo.R`, that can be run without
  reading the lower-level function files.
- A user-facing function, `run_sim2_demo()`, that wraps the full optimized
  local pipeline.
- A simplified flat `R/` folder with broad files for setup, data generation,
  model fitting, AIPW estimation, and reporting.
- `demo` sample scale with values `120` and `300`.
- CLI flags for local CPU/CV controls such as `--total-cores`,
  `--reserve-cores`, `--fit-workers`, `--cv-workers-per-fit`,
  `--max-concurrent-hal-fits`, `--max-active-mc`, `--cv-folds`, and
  `--cv-nlambda`.
- Lightweight local report generation plus the three standard Simulation 2
  plots: estimates, variance calibration, and confidence interval coverage.
- A one-worker CV fast path that avoids `makeCluster(1)`.
- A cross-run global queue that flattens complete HAL fits across every selected
  sample-size-combination-by-MC run while retaining a per-run fit cap.
- Persistent `callr` HAL slots. With one CV worker a slot fits directly; with
  more than one it owns a private reusable CV PSOCK pool that never crosses a
  process boundary.
- Immediate Naive/Proposed AIPW assembly, checkpointing, and memory release as
  soon as all four nuisance fits for one run succeed.
- A bounded active-run window that limits how many prepared, unfinished run
  contexts and their large data/fit objects can be resident at once.
- Content-addressed, validated, atomic HAL fit caches.
- Run-level checkpoints, resume, and per-run atomic CSV shards.
- A scheduler event trace plus separate preparation, first-queue-wait, fit-span,
  and end-to-end invocation-wall timing fields.
- Local cache and Brian-style output directories inside `aipw_local_demo_2/`.

The adopted automatic scheduling notation is:

```text
U = usable CPU units after reservation
c = CV workers per complete HAL fit (default 1)
F = simultaneous HAL fits allowed per run (default 2)
S = floor(U / c) global HAL-slot CPU ceiling
A = min(R, ceiling(1.2 * S / F)) active-run memory window
```

Here `R` counts sample-size-combination-by-MC runs. Effective `S` is clipped
when the workload, an explicit slot limit, or a manually reduced active window
cannot feed the CPU-derived ceiling. The direct CLI without a profile reserves
one CPU unit when possible; the named profiles explicitly reserve zero.

These queue, worker, and persistence choices do not change the DGP, fold IDs,
four nuisance-fit roles, or AIPW equations. They are distinct from the current
three-fold and `num_knots = c(5, 3)` tuning reductions documented above, which
do change the fitted nuisance models.

## Intended Use

Use this folder for teaching, code review, laptop checks, and future GitHub
tutorial material. The active-run bound protects coordinator memory but cannot
make an arbitrarily large sample fit in limited RAM: each persistent slot and
any private CV workers are separate R processes, while each active run retains
analysis data, four training inputs, and completed HAL objects until immediate
finalization. Lower `--max-active-mc` and then
`--max-concurrent-hal-fits` on constrained machines. Use the existing
Slurm-oriented project folders for full Monte Carlo production runs.
