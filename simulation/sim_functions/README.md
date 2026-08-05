# Simulation function location

The canonical `run_sim2_demo()` and reporting implementations now live in
`R/sim2_demo.R` and `R/sim2_demo_report.R` so they are included when the
`aipwLocalDemo` package is installed. The repository CLI continues to load
those files through `R/source_all.R`.
