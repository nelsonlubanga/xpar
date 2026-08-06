# Reproducible simulation study backing vignette("02-simulation-validation").
# Four sweeps around a shared baseline condition (divergence = 2000
# generations, 10% missingness, 10 candidates/species, error = 0), each
# varying one axis at a time: genotyping error rate, divergence, missingness,
# candidate-pool size. Re-run with `source("data-raw/simulation_study.R")`
# from the package root; takes several minutes (AlphaSimR + block bootstrap
# per trial).

devtools::load_all(".", quiet = TRUE)

# 1. Error rate x divergence (broad) ------------------------------------
sweep_error_broad <- benchmark_grid(
  divergence = c(20, 2000), error_rate = c(0, 0.10, 0.25, 0.40),
  missing_rate = 0.10, n_per_species = 10, n_chr = 10, seg_sites = 200,
  n_boot = 300, n_rep = 5, seed = 200
)

# 2. Error rate, fine-grained around the breakdown point (high divergence) -
sweep_error_fine <- benchmark_grid(
  divergence = 2000, error_rate = c(0.25, 0.30, 0.33, 0.35, 0.38, 0.40),
  missing_rate = 0.10, n_per_species = 10, n_chr = 10, seg_sites = 200,
  n_boot = 300, n_rep = 10, seed = 300
)

# 3. Missingness, at a fixed moderate error rate -------------------------
sweep_missing <- benchmark_grid(
  divergence = 2000, error_rate = 0.20,
  missing_rate = c(0, 0.10, 0.20, 0.30, 0.40),
  n_per_species = 10, n_chr = 10, seg_sites = 200,
  n_boot = 300, n_rep = 8, seed = 400
)

# 4. Candidate-pool size, at a fixed moderate error rate -----------------
sweep_pool_size <- benchmark_grid(
  divergence = 2000, error_rate = 0.20, missing_rate = 0.10,
  n_per_species = c(3, 5, 10, 20), n_chr = 10, seg_sites = 200,
  n_boot = 300, n_rep = 8, seed = 500
)

sim_study <- list(
  error_broad = sweep_error_broad,
  error_fine = sweep_error_fine,
  missingness = sweep_missing,
  pool_size = sweep_pool_size
)

usethis::use_data(sim_study, overwrite = TRUE)
