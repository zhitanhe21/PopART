# Function Map - AIPW Local Demo 2

The user-facing scripts stay small; computational, persistence, and reporting
helpers live in the package `R/` files.

## User-facing workflow

| Location | Role |
|---|---|
| `simulation/sim_scripts/ss2_demo.R` | Main RStudio/CLI entry point; selects local defaults, runs the simulation, and optionally writes the report and figures. |
| `R/profiles.R` | Defines the optional `smoke`, `quick`, `hour`, and `formal` presets plus `run_sim2_demo_profile()`. |
| `R/sim2_demo.R` | Defines `run_sim2_demo()`, including pool lifecycle, resume, and aggregate output. |
| `simulation/sim_analysis/sim2_demo_analysis.R` | Independently reads existing CSVs and creates the report and figures. |
| `R/sim2_demo_report.R` | Defines the report and three standard Simulation 2 figures. |

## Implementation files

| File | Role |
|---|---|
| `R/source_all.R` | Loads functions in dependency order. |
| `R/profiles.R` | Resolves named workload defaults; explicit arguments take precedence. |
| `R/setup_helpers.R` | Package loading, CLI helpers, sample grids, and local CPU planning. |
| `R/persistence.R` | Content hashing, validated fit caches, immutable run checkpoints/shards, and atomic file promotion. |
| `R/sim2_data.R` | Preserved Simulation 2 data generation. |
| `R/model_fitting.R` | HAL fitting, CV fast path, timing, reusable PSOCK pools, and load-balanced dispatch. |
| `R/aipw_estimators.R` | Four nuisance-model jobs plus naive/proposed AIPW estimation. |
| `R/report_helpers.R` | Lightweight Markdown table formatting. |

## Execution flow

```text
simulation/sim_scripts/ss2_demo.R
  -> optional sim2_profile_defaults()
  -> run_sim2_demo()
       -> fitcv_sample_grid()
       -> configure_fitcv_parallel()
       -> restore a matching successful checkpoint, if available
       -> create one reusable PSOCK pool, if computation is needed
       -> sim2_fun_fitcvreuse_local()
            -> simulate_sim2_data()
            -> make_fitcvreuse_jobs()
            -> run_hal_jobs_local()
            -> compute_naive_eta()
            -> compute_proposed_eta()
       -> atomically promote immutable checkpoint and per-run shards
       -> recoverably rebuild aggregate CSVs and sd0.csv
  -> write_sim2_demo_report()
  -> write_sim2_demo_figures()

simulation/sim_analysis/sim2_demo_analysis.R
  -> write_sim2_demo_report()
  -> write_sim2_demo_figures()
       -> *_estimates.png
       -> *_variance.png
       -> *_confidence.png

Programmatic profile use:

```text
run_sim2_demo_profile()
  -> resolve_sim2_profile_args()
  -> run_sim2_demo()
```
```
