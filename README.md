# aipwLocalDemo

This directory is an independent version-2 copy of `aipw_local_demo`. It keeps
the Simulation 2 DGP, nuisance-job structure, and AIPW formulas. The current
experimental laptop defaults use three-fold CV and `num_knots = c(5, 3)`, so
the fitted HAL models are not numerically identical to the original settings.
The local execution layer adds:

- a global queue of complete HAL fits spanning every selected sample-size
  combination and Monte Carlo replicate;
- persistent `callr` HAL slots that take the next eligible fit as soon as they
  become free;
- a true serial fast path when `cv_workers_per_fit = 1`, plus one private,
  reusable CV PSOCK pool per `callr` slot when more CV workers are requested;
- content-addressed, validated, atomic HAL fit caches;
- successful-run checkpoints and automatic resume;
- per-run atomic CSV shards plus final compatible aggregate CSV files.

The original `aipw_local_demo` directory is not modified.

## Install as an R package

From a local checkout:

```bash
R CMD INSTALL aipw_local_demo_2
```

Then use the package from any working directory:

```r
library(aipwLocalDemo)

plan <- run_sim2_demo(
  run_id = "local_plan",
  max_runs = 1,
  total_cores = 4,
  reserve_cores = 1,
  dry_run = TRUE
)
```

The GitHub repository is initially private because the licensing terms for the
upstream research code still need to be confirmed. See `LICENSE` and
`docs/provenance.md` before redistributing the source.

## Quick start

From RStudio:

1. Open `aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R`.
2. Source the complete script.

From a terminal:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R
```

The global scheduler requires the R package `callr`, in addition to the HAL and
reporting packages used by the original local demo.

The default CLI run remains intentionally small: one 120/120 sample-size pair,
three CV folds, `num_knots = c(5, 3)`, and five lambda values. The programmatic
API keeps the research default of 50 lambda values unless it is overridden.

The automatic global plan uses the following quantities:

```text
R = number of selected sample-size-combination x MC runs
U = detected/overridden CPU units minus reserved units
c = CV workers owned by one HAL fit (default 1)
F = maximum simultaneous HAL fits from one run (default 2)
S = floor(U / c) global HAL-slot CPU ceiling
A = min(R, ceiling(1.2 * S / F)) prepared-run memory window
```

The reported effective slot count is additionally clipped when there are fewer
than `S` eligible fits, when `--max-concurrent-hal-fits` is lower, or when a
manually reduced active window can feed only `A * F` fits. Without a profile,
the CLI reserves one CPU unit when possible; named profiles explicitly reserve
zero. With four detected units and the one-run tiny default, `U = 3`, `c = 1`,
`F = 2`, `A = 1`, and only two of the three CPU-derived slots are needed. No
inner CV process is created when `c = 1`.

## Validation profiles

The historical 5.44-hour benchmark is the publication-style `4 sample
combinations x 20 MC` workload. It was measured with the previous per-run
queue, four-fold CV, and `num_knots = c(25, 10)` configuration and is not a
timing estimate for the current global scheduler and experimental defaults.
Four opt-in profiles make the intended workload explicit:

| Profile | Workload | Current CV controls | Intended use | Historical 4-fold / c(25,10) timing evidence |
|---|---:|---:|---|---:|
| `smoke` | one 120/120 run | 3 folds, 50 lambdas | Check that the code and packages run. | 2.65 s recorded cold run; 4.074 s invocation wall on the benchmark laptop. |
| `quick` | four small-grid runs, 1 MC | 3 folds, 50 lambdas | Validate every current sample-size path. | 12.80 min for the first formal MC; all 20 observed MCs were 10.98--25.80 min. |
| `hour` | eight small-grid runs, 2 MC | 3 folds, 50 lambdas | Repeat the current replicate algorithm with two seeds. | First two MCs: 30.17 min; every pair among the 20 observed MCs was at most 50.83 min. |
| `formal` | 80 small-grid runs, 20 MC | 3 folds, 50 lambdas | Full statistical experiment under the current tuning configuration. | 5.4423 h actual. |

Run a profile with:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --profile quick \
  --no-analysis \
  --no-plots
```

Use `--profile smoke` for the shortest check, `--profile hour` for two MC
replicates, or `--profile formal` for the complete benchmark. Explicit CLI
arguments override profile defaults. For example, one `(500, 500)` formal-
configuration replicate is:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --profile quick \
  --max-runs 1 \
  --no-analysis \
  --no-plots
```

Under the historical four-fold, `c(25, 10)` configuration, that combination
averaged 43.86 seconds in the formal run. A single
`(5000, 5000)` replicate can be selected with
`--n-trial 5000 --n-aux 5000`; it averaged 7.11 minutes and had an observed
maximum of 11.61 minutes.

`quick` and `hour` reduce only the number of Monte Carlo replicates. Their
retained DGP, folds, lambda grid, HAL fits, and AIPW calculations match
`formal`. They are appropriate for execution and numerical regression checks,
but one or two MC replicates are not enough for stable bias, variance, or
coverage conclusions.

Profiles are optional presets. Running the CLI without `--profile` preserves
the tiny default described above, including `reserve_cores = 1`. Profiles use
`reserve_cores = 0`; `smoke` uses `F = 1`, while `quick`, `hour`, and `formal`
use the default per-run cap `F = 2`.

## Recommended laptop command

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --scale demo \
  --max-runs 1 \
  --total-cores 4 \
  --reserve-cores 1 \
  --fit-workers 2 \
  --cv-workers-per-fit 1 \
  --cv-nlambda 5
```

Worker counts are optional; omitting them uses the version-2 automatic plan.
The detector prefers physical cores where the operating system exposes them,
but may receive logical processors or a scheduler allocation instead. Treat
the printed `total_cores` value as a best-effort CPU budget, not a guarantee
about hardware topology; use `--total-cores` to set a known-safe budget.

## Monte Carlo grid

`--mc-reps N` runs every selected sample-size combination `N` times. As in
Brian's original Simulation 2 grid, replicate `r` uses `base_seed + r` for all
selected sample sizes. `--max-runs` selects sample-size combinations before
this expansion. The resulting runs share one global HAL queue, so a slot freed
by a fit for one sample size or MC replicate may immediately start an eligible
fit from another one.

For the complete `(500, 5000)` four-combination grid with 20 MC replicates:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --scale small \
  --max-runs 4 \
  --mc-reps 20 \
  --total-cores 4 \
  --reserve-cores 0 \
  --fit-workers 2 \
  --cv-workers-per-fit 1 \
  --cv-folds 3 \
  --cv-nlambda 50 \
  --no-cache \
  --no-analysis \
  --run-id small_mc20_globalS4_F2_cv1_folds3_knots5_3
```

This is an 80-run compute benchmark. With four usable CPU units and `c = 1`,
the automatic plan has `S = 4`, `F = 2`, and
`A = ceiling(1.2 * 4 / 2) = 3`. Increase `reserve_cores`, reduce
`--max-concurrent-hal-fits`, or reduce `--max-active-mc` for a more responsive
or memory-constrained laptop run. Checkpoints remain enabled, so rerunning the
same command resumes completed rows.

## Global HAL scheduler

Each run represents one sample-size combination and one MC replicate. After a
run is admitted, the coordinator simulates its data and constructs the four
complete nuisance fits:

```text
proposed Q0   proposed Q1   shared outcome   naive censor
       \          |              |              /
        +--------- global eligible-fit queue --+
                          |
              S persistent callr slots
```

At most `F` fits from one run may execute simultaneously. When any persistent
slot finishes, it takes the next eligible fit globally; it does not wait for the
other fits in the same run. As soon as all four fits for a run succeed, the main
R process immediately computes its Naive and Proposed AIPW rows, writes its
atomic checkpoint and shards, releases that run's large data and fit objects,
and admits another run if work remains. It does not wait for the other active
runs or for every sample size sharing the same MC seed.

With `c = 1`, a `callr` slot directly performs the complete HAL fit and the
`cv.glmnet` folds run serially inside it. With `c > 1`, each persistent slot
owns a private reusable CV PSOCK pool of size `c`; socket-cluster objects never
move between processes. Increasing `c` therefore reduces the number of global
fit slots and increases process and memory overhead.

`A` is a memory valve, not another CPU multiplier. An active context can retain
the combined analysis data, four training matrices and responses, and completed
HAL objects until its fourth fit finishes. On a memory-limited PC, lower
`--max-active-mc` first; if necessary also lower
`--max-concurrent-hal-fits`. Large sample grids can still exceed laptop memory
even when the CPU equations are satisfied.

## Timing interpretation

- `*_run_summary.csv: invocation_wall_elapsed_sec` is the end-to-end waiting
  time for the invocation and is the correct primary scheduler/runtime measure.
- `*_timing.csv: prepare_elapsed_sec` measures admission through preparation.
- `queue_wait_sec` measures preparation through the first dispatched fit; it is
  not the sum of all later waits between the four fits.
- `fit_span_sec` measures first dispatch through the last completed fit and can
  include staggered queue gaps.
- The run-level `elapsed_sec` covers admission through AIPW assembly after all
  four fits complete. The subsequent checkpoint/shard writes remain part of
  invocation wall time but are not currently split into their own per-run
  timing column.
- `recorded_run_elapsed_sec` is a sum of overlapping per-run latencies under the
  global scheduler. It is work/latency accounting, not wall time, and must not
  be used as the global speed comparison.

The `*_scheduler_trace.csv` event log and the observed maxima in the run summary
show whether the global, per-run, and active-context caps were respected. A
fully resumed invocation writes `resume` events and correctly reports observed
concurrency maxima of zero because it starts no scheduler slots.

## Cache and resume

Run checkpoints are enabled by default. If a complete run has a valid matching
checkpoint, restarting the same command skips that run. Use `--no-resume` when
you intentionally want to recompute it.

HAL fit caching is separate. It is opt-in in the direct CLI because full fitted
objects can consume substantial disk space:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R --use-cache
```

The programmatic `run_sim2_demo()` API retains the inherited
`use_cache = TRUE` default; pass `use_cache = FALSE` if disk use matters.

The cache key includes the actual training `X`, `Y`, weights, fold IDs, HAL/CV
controls, and relevant package/R versions. CPU counts, worker counts, and
backend are excluded because they change scheduling rather than the statistical
fit. Invalid cache files are moved aside and recalculated.

Two recovery levels are therefore available:

```text
fit cache       -> recover individual completed HAL fits
run checkpoint  -> skip a complete successful sample-grid run
```

## Atomic outputs

Each successful run first writes same-directory temporary files and promotes
them to a checkpoint and three per-run CSV shards:

```text
<out_dir>/checkpoints/<prefix>/
  run_0001_<signature>_<attempt-id>.rds
  csv/
    run_0001_<signature>_<attempt-id>_timing.csv
    run_0001_<signature>_<attempt-id>_fit_timing.csv
    run_0001_<signature>_<attempt-id>_results.csv
```

After all requested runs finish or resume, the program rebuilds the compatible
aggregate files through a same-directory temporary file and a recoverable
replacement:

```text
<prefix>_parallel_plan.csv
<prefix>_timing.csv
<prefix>_fit_timing.csv
<prefix>_results.csv
<prefix>_run_summary.csv
<prefix>_scheduler_trace.csv  # scheduler events, or resume-only events
sd0.csv
```

The checkpoints and shards are the recovery source; the aggregate CSV files are
derived reporting files.

`--no-resume` preserves earlier attempts and creates a new attempt-id, so these
immutable recovery files can accumulate. A failed attempt always retains its
checkpoint and timing shard; fit/result shards exist only when those tables
were produced.

## Useful options

| Option | Meaning |
|---|---|
| `--profile smoke\|quick\|hour\|formal` | Apply a named workload preset; explicit options override it. |
| `--fit-workers N` | Per-run HAL-fit cap `F`; the automatic default is 2. |
| `--cv-workers-per-fit N` | CV worker processes owned by each fit. |
| `--max-concurrent-hal-fits N` | Cap the effective global slots; it cannot exceed the CPU ceiling `S`. |
| `--max-active-mc N` | Bound prepared, unfinished run contexts; lower it to reduce coordinator memory. |
| `--active-mc-headroom X` | Set the multiplier in automatic `A`; default `1.2`. |
| `--backend auto\|psock\|fork` | Compatibility field; global execution uses `callr`, and `fork` remains invalid on Windows. |
| `--mc-reps N` | MC replicates per selected sample-size combination. |
| `--max-runs N` | Sample-size combinations selected before MC expansion. |
| `--use-cache` | Enable validated HAL fit caching. |
| `--no-cache` | Disable HAL fit caching. |
| `--no-resume` | Ignore matching successful run checkpoints and recompute. |
| `--checkpoint-dir PATH` | Override the run checkpoint directory. |
| `--no-analysis` | Skip the Markdown report and plots. |
| `--no-plots` | Write the report but skip plots. |
| `--dry-run` | Print the plan without creating `callr` slots or fitting HAL. |
| `--help` | Show all primary CLI options without starting a run. |

For a manually selected sequential-fit/two-CV-worker layout:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --fit-workers 1 \
  --cv-workers-per-fit 2
```

In this layout, `F = 1`. Each active global `callr` slot owns a private
two-worker CV PSOCK pool and reuses that private pool for every fit assigned to
the slot.

The original Simulation 2 data-generating intercepts are fixed for
`p_resp = 0.5` and `p_cens = 0.3`. Version 2 rejects other values rather than
writing misleading metadata; supporting other targets would require deriving
new intercepts and would be an algorithm change.

Do not run two processes against the same `out_dir`. Immutable attempt files
have unique names, but aggregate CSVs (including `sd0.csv`) are intentionally
single-writer files and are not protected by a cross-process lock. Give each
concurrent invocation its own output directory.

## Main files

- `R/profiles.R`: named workload presets and the programmatic profile wrapper.
- `R/setup_helpers.R`: local CPU planning and command-line parsing.
- `R/persistence.R`: hashing, cache validation, checkpoints, and atomic files.
- `R/model_fitting.R`: one complete HAL fit, CV fast path, and fit caching.
- `R/global_scheduler.R`: bounded cross-run queue and persistent `callr` slots.
- `R/aipw_estimators.R`: four nuisance tasks and the two AIPW estimators.
- `R/sim2_demo.R`: global-run admission, immediate
  four-fit finalization, resume, and aggregate output orchestration.
- `R/sim2_demo_report.R`: reports and the three standard diagnostic figures.
- `tests/test_global_scheduler.R`: scheduler limits, cross-run refill, and
  failure isolation with lightweight tasks.
- `tests/test_persistence_helpers.R`: lightweight planning/persistence checks.
- `tests/test_profiles.R`: profile resolution, dry-run grids, and signature
  equivalence.

## Validation

Run the lightweight checks with:

```bash
Rscript aipw_local_demo_2/tests/test_persistence_helpers.R
Rscript aipw_local_demo_2/tests/test_profiles.R
Rscript aipw_local_demo_2/tests/test_global_scheduler.R
```

Run the package build checks with:

```bash
R CMD build aipw_local_demo_2
R CMD check --no-manual aipwLocalDemo_0.1.0.tar.gz
```

For numerical regression, use the same seed, folds, sample size, and HAL
controls in both demo versions and compare the estimator columns. Worker
scheduling should not change them, apart from possible platform-level floating
point differences.
