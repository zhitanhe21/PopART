#!/usr/bin/env bash
set -euo pipefail

# Run from the HAL project root. This follows the Brian-style layout.
Rscript aipw_local_demo_2/simulation/sim_scripts/ss2_demo.R \
  --scale demo \
  --max-runs 1 \
  --run-id demo_run \
  --reserve-cores 2 \
  --cv-workers-per-fit 1 \
  --cv-nlambda 5
