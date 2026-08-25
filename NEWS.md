# popart 0.2.0

- Added independent learner selection for the outcome, censoring, Q0, and Q1
  nuisance models.
- Added XGBoost nuisance fitting with explicit boosting controls, reproducible
  seeds, and one thread per complete fit.
- Preserved HAL as the default learner for all nuisance models.

# popart 0.1.0

- Added default three-fold outer cross-fitting with pooled out-of-fold nuisance
  predictions.
- Added CPU-aware complete-fit scheduling, HAL internal cross-validation,
  diagnostics, documentation, and repository simulation tools.
