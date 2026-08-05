# Local Demo 2 Workflow

This guide describes the version-2 local AIPW/HAL workflow. The DGP, estimator
formulas, and four nuisance-job roles are retained. The current experimental
laptop configuration also changes HAL tuning to three-fold CV and
`num_knots = c(5, 3)`, in addition to the scheduling and persistence changes.

## What to run

From RStudio:

1. Open `aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R`.
2. Source the complete script.

From the workspace root:

```bash
bash aipw_local_demo_2/examples/run_demo.sh
```

or:

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R
```

The independent analysis entry point is
`aipw_local_demo_2/simulation/sim_analysis/sim2_demo_analysis.R`.

## Choose a workload, not a five-hour default

The direct CLI still has its original tiny default. Named profiles are
optional argument presets:

| Profile | Sample grid x MC | Folds / lambdas | Typical purpose |
|---|---:|---:|---|
| `smoke` | one demo run | 3 / 50 | Package and code-path check. |
| `quick` | four small-grid runs x 1 MC | 3 / 50 | Complete current-configuration path check. |
| `hour` | four small-grid runs x 2 MC | 3 / 50 | Two-replicate validation. |
| `formal` | four small-grid runs x 20 MC | 3 / 50 | Full 80-run statistical experiment under the current tuning. |

```bash
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --profile quick --no-analysis --no-plots
```

Historical timing evidence for one complete four-combination MC ranged from
10.98 to 25.80 minutes over 20 seeds. Those measurements used four folds and
`num_knots = c(25, 10)` under the previous per-run queue; they are retained as
provenance, not as predictions for the current global scheduler and tuning.
Other applications, cooling, CPU generation, and available memory can also
change wall time.

`quick` and `hour` change MC precision, not the calculation within a retained
replicate. Do not use one or two replicates to draw stable Monte Carlo bias or
coverage conclusions. Explicit options always override a profile. Omitting
`--profile` keeps the tiny CLI workload and reserves one CPU unit when
possible; profiles explicitly reserve zero CPU units.

## Default plan

| Setting | CLI default |
|---|---|
| Sample scale | `demo` |
| First sample-size pair | `120 / 120` |
| Maximum grid rows | `1` |
| MC replicates per selected size | `1` |
| CPU budget | Best-effort auto-detected |
| Reserved system CPU units | `1` when possible without a profile; profiles use `0` |
| Per-run complete-fit cap `F` | `2` (`smoke` uses `1`) |
| CV workers per fit | `1` |
| Global HAL-slot CPU ceiling `S` | `floor(U / c)` |
| Prepared-run window `A` | `min(R, ceiling(1.2 * S / F))` |
| CV folds | `3` |
| HAL knots by degree | `c(5, 3)` |
| Lambda values | `5` |
| Run resume | Enabled |
| Full HAL fit cache | Disabled unless `--use-cache` is supplied |

Here, `R` is the number of selected sample-size-combination-by-MC runs, `U` is
the usable CPU budget after reservation, and `c` is the number of CV workers
owned by one complete HAL fit. `S = floor(U / c)` is the core-derived ceiling;
the effective value printed as `max_concurrent_hal_fits` is also clipped by
available work (`R * F`), explicit overrides, and a manually reduced `A * F`.
With a detected budget of four CPU units and no profile, `U = 3` and the
one-run tiny workload needs two slots because its per-run cap is `F = 2`.
`c = 1` runs CV in the persistent fit process and creates no inner socket
worker.

Detection prefers physical cores when the operating system exposes them, but
can fall back to logical processors or a scheduler allocation. The printed
`total_cores` value is therefore a planning budget, not a guarantee about
physical topology. Set `--total-cores` explicitly when you know the safe
budget for the computer.

## Cross-run global queue

The global scheduler flattens eligible complete HAL fits across sample sizes
and Monte Carlo replicates:

```text
active run 1: Q0 Q1 shared-outcome censor --+
active run 2: Q0 Q1 shared-outcome censor --+--> global eligible queue
active run 3: Q0 Q1 shared-outcome censor --+             |
                                                       S persistent
                                                       callr slots
```

At most `F` jobs from one run can be active. The first free slot takes the next
eligible fit from any admitted run, so it can cross both sample-size and MC
boundaries rather than waiting behind a slow fit in the same run. The `callr`
sessions persist for the invocation. With `c > 1`, each session owns its own
private reusable CV PSOCK pool; a socket-cluster object is never transferred
between processes. The package `callr` is therefore required.

For Monte Carlo work, `--mc-reps N` expands each selected sample-size
combination into `N` runs. Replicate `r` uses `base_seed + r` across all sample
sizes, matching Brian's original `expand.grid(..., sim.id=...)` seed layout.
For example, `--scale small --max-runs 4 --mc-reps 20` creates 80 runs.

`A` limits prepared, unfinished run contexts and is primarily a memory valve.
Each context can retain its analysis data, four training inputs, and completed
HAL objects. Once all four fits for one run succeed, the main R process
immediately computes both AIPW versions, writes the atomic run checkpoint and
shards, removes that context, and admits another run. It does not wait for the
other active runs or for all sample sizes with the same MC seed.

For a RAM-constrained laptop, lower `--max-active-mc`; if needed, also lower
`--max-concurrent-hal-fits`. Increasing `c` creates more worker processes per
slot and reduces `S`, while large samples and HAL dictionaries can remain too
large for a laptop regardless of the CPU plan.

## Persistence flow

```text
build sample-size x MC run identities
  |
  +-- restore matching successful checkpoints before preparing data
  |
  +-- admit at most A missing runs and enqueue their four HAL fits
        |
        +-- at most S fits globally and F fits per run
        +-- optional validated fit cache per HAL job
        +-- fourth successful fit -> immediate AIPW assembly
              |
              +-- atomically promote run checkpoint and CSV shards
              +-- release context and admit the next run
                    |
                    +-- rebuild aggregate CSV files after all runs
```

Fit-cache identity is calculated from the actual training data, weights, folds,
HAL/CV controls, and relevant versions. Worker layout is excluded, allowing a
fit created under `2 x 1` to be reused under `1 x 2`.

Run-checkpoint identity also includes the seed, requested sample sizes,
response/censoring settings, computational-source fingerprint, and statistical
controls. Only successful checkpoints containing two estimator rows and four
fit-timing rows are eligible for resume. Editing tests, analysis scripts, or
report formatting does not invalidate a compatible computational checkpoint.

Use:

```bash
--use-cache      # keep and reuse complete HAL fit objects
--no-resume      # intentionally recompute a completed run
```

## Outputs

Aggregate report-compatible files remain under:

```text
aipw_local_demo_2/simulation/sim_data/sim2/
aipw_local_demo_2/simulation/sim_figures/sim2/
```

| File pattern | Meaning |
|---|---|
| `sd0.csv` | Brian-style estimator output. |
| `*_results.csv` | Aggregated estimator results. |
| `*_timing.csv` | Aggregated run timing and status. |
| `*_fit_timing.csv` | Aggregated four-fit diagnostics. |
| `*_parallel_plan.csv` | Selected worker plan. |
| `*_run_summary.csv` | End-to-end invocation wall time and row counts. |
| `*_scheduler_trace.csv` | Admission, dispatch, completion, and cap-observation events; a fully restored invocation records resume-only events. |
| `checkpoints/<prefix>/*_<attempt-id>.rds` | Immutable complete-run recovery checkpoints. |
| `checkpoints/<prefix>/csv/*_<attempt-id>_*.csv` | Immutable per-run output shards. |

The checkpoints and shards are the durable source. Aggregate CSV files are
rebuilt after all requested runs finish or resume using same-directory
temporary files and a recoverable replacement.

### Timing fields

- `*_run_summary.csv: invocation_wall_elapsed_sec` is the full user-visible
  waiting time and is the correct primary measure for global-scheduler speed.
- In `*_timing.csv`, `prepare_elapsed_sec` is admission through preparation,
  `queue_wait_sec` is preparation through the first dispatch, and
  `fit_span_sec` is first dispatch through the fourth completed fit. The latter
  can include staggered waits between fits.
- Run-level `elapsed_sec` ends after AIPW assembly once all four fits complete.
  The subsequent checkpoint/shard writes are included in invocation wall time
  but are not currently a separate per-run field.
- `recorded_run_elapsed_sec` sums overlapping run latencies and is not elapsed
  wall time under a global queue. Do not use it for scheduler speedups.
- On a fit-cache hit, fit-level `elapsed_sec` retains the original fit cost for
  provenance, while `effective_elapsed_sec` records the current lookup cost.

The scheduler trace and `max_observed_global_fits`,
`max_observed_per_run_fits`, and `max_observed_active_runs` in the summary can
be used to audit the three concurrency limits. All three observed maxima are
zero for a fully resumed invocation because it starts no worker slots.

## Programmatic use

```r
project_dir <- normalizePath("aipw_local_demo_2")
source(file.path(project_dir, "R", "source_all.R"))
source_aipw_local_demo(project_dir)

run_sim2_demo(
  project_dir = project_dir,
  run_id = "demo_run",
  max_runs = 1,
  mc_reps = 1,
  total_cores = 4,
  reserve_cores = 1,
  fit_workers = 2,
  cv_workers_per_fit = 1,
  cv_nlambda = 5,
  use_cache = FALSE,
  resume = TRUE,
  write_sd_file = TRUE
)
```

The same presets are available programmatically:

```r
run_sim2_demo_profile(
  project_dir = project_dir,
  profile = "quick",
  run_id = "quick_validation",
  use_cache = FALSE
)
```

Profile names are convenience metadata only. Resolved statistical arguments,
not the profile name, determine cache and checkpoint compatibility.

The programmatic API defaults to 50 lambda values, preserving the research
setting. The direct CLI defaults to 5 to keep an initial laptop run short.

Simulation 2 also preserves its original fixed response/censoring targets:
`p_resp = 0.5` and `p_cens = 0.3`. Other values would require deriving new DGP
intercepts, so version 2 rejects them instead of treating them as tuning
options.

The persistence layer assumes one writer per output directory. For concurrent
invocations, use a separate `out_dir` for every process; aggregate files and
`sd0.csv` do not use a cross-process lock.

## Main implementation

| Function | File | Role |
|---|---|---|
| `run_sim2_demo_profile()` | `R/profiles.R` | Resolve a named workload preset and forward explicit overrides. |
| `run_sim2_demo()` | `R/sim2_demo.R` | Global-run admission, checkpoint restore, immediate per-run finalization, and aggregate outputs. |
| `run_global_hal_scheduler()` | `R/global_scheduler.R` | Maintain the bounded active-run window and dispatch eligible fits to persistent `callr` slots. |
| `fit_hal_timed()` | `R/model_fitting.R` | HAL fit, CV fast path, fit cache, and timing. |
| `run_one_hal_job()` | `R/model_fitting.R` | Execute one complete nuisance fit inside a global slot. |
| `make_fit_cache_identity()` | `R/persistence.R` | Statistical cache identity. |
| `atomic_save_rds()` | `R/persistence.R` | Same-directory temporary write and promotion. |
| `find_successful_run_checkpoint()` | `R/persistence.R` | Successful immutable-checkpoint lookup. |
| `sim2_fun_fitcvreuse_local()` | `R/aipw_estimators.R` | One complete simulation run. |
