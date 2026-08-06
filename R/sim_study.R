#' Simulation study: recovery power under controlled conditions
#'
#' Pre-computed results from [benchmark_grid()] sweeps used in
#' `vignette("02-simulation-validation")`, characterizing how often `xpar`
#' recovers the true parent pair as conditions worsen. Regenerate with
#' `source("data-raw/simulation_study.R")` (several minutes; requires
#' \pkg{AlphaSimR}).
#'
#' All four sweeps share a baseline condition (divergence 2000 generations,
#' 10 candidates per species, 10 linkage groups, 200 segregating sites per
#' chromosome) and vary one axis at a time.
#'
#' @format A list of four data frames, each with one row per
#'   condition x replicate as returned by [benchmark_grid()] (columns:
#'   `divergence`, `error_rate`, `missing_rate`, `n_per_species`, `n_pairs`,
#'   `true_pair_rate`, `true_pair_rank`, `true_pair_is_top`,
#'   `true_pair_win_prob`, `rep`):
#'   \describe{
#'     \item{error_broad}{Genotyping error rate (0-40%) crossed with low
#'       (20 gen) and high (2000 gen) divergence.}
#'     \item{error_fine}{Genotyping error rate, finely spaced (25-40%) at
#'       high divergence only, to localize the recovery breakdown point.}
#'     \item{missingness}{Missing-call rate (0-40%) at a fixed 20% error
#'       rate and high divergence.}
#'     \item{pool_size}{Candidate-pool size per species (3-20) at a fixed
#'       20% error rate, 10% missingness, and high divergence.}
#'   }
#' @source `data-raw/simulation_study.R`
"sim_study"
