#' Run one simulated recovery trial
#'
#' Simulates a two-species candidate pool and known cross
#' ([simulate_hybrid_pedigree()]), overlays noise ([inject_noise()]), runs
#' [block_bootstrap()] over every candidate pair, and scores whether the
#' method recovered the true parent pair. This is the unit of replication
#' for [benchmark_grid()].
#'
#' @param n_per_species,divergence,n_chr,seg_sites Passed to
#'   [simulate_hybrid_pedigree()].
#' @param error_rate,missing_rate Passed to [inject_noise()].
#' @param n_boot Passed to [block_bootstrap()]. Default `500` (lower than
#'   the [block_bootstrap()] default, since [benchmark_grid()] calls this
#'   many times).
#' @param seed Optional integer seed for reproducibility. Governs both the
#'   simulation and the noise injection.
#'
#' @return A one-row data frame: `divergence`, `error_rate`, `missing_rate`,
#'   `n_per_species`, `n_pairs`, `true_pair_rate`, `true_pair_rank`,
#'   `true_pair_is_top` (logical — was the true pair the single best-ranked
#'   pair), `true_pair_win_prob`.
#' @export
run_trial <- function(n_per_species = 5, divergence = 1000, n_chr = 10,
                       seg_sites = 200, error_rate = 0, missing_rate = 0,
                       n_boot = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  sim <- simulate_hybrid_pedigree(
    n_per_species = n_per_species, divergence = divergence,
    n_chr = n_chr, seg_sites = seg_sites
  )
  noisy_geno <- inject_noise(sim$geno, error_rate = error_rate, missing_rate = missing_rate)
  d <- xpar_data(noisy_geno, sim$scaffold)

  boot <- block_bootstrap(
    d, candidates = sim$candidates, offspring = sim$offspring,
    self = FALSE, n_boot = n_boot
  )
  ranked <- rank_pairs(boot)
  ranked$rank <- seq_len(nrow(ranked))

  is_true <- (ranked$parent1 == sim$true_parent1 & ranked$parent2 == sim$true_parent2) |
    (ranked$parent1 == sim$true_parent2 & ranked$parent2 == sim$true_parent1)
  true_row <- ranked[is_true, ]

  data.frame(
    divergence = divergence,
    error_rate = error_rate,
    missing_rate = missing_rate,
    n_per_species = n_per_species,
    n_pairs = nrow(ranked),
    true_pair_rate = true_row$rate,
    true_pair_rank = true_row$rank,
    true_pair_is_top = true_row$rank == 1,
    true_pair_win_prob = true_row$win_prob
  )
}

#' Benchmark recovery power over a grid of simulation conditions
#'
#' Runs [run_trial()] over every combination of the supplied parameter
#' vectors, with `n_rep` independent replicates per combination, and returns
#' the stacked results — the basic input for a power / false-assignment-rate
#' summary (e.g. `mean(true_pair_is_top)` per condition).
#'
#' @param divergence,error_rate,missing_rate,n_per_species Vectors of values
#'   to cross (via [expand.grid()]) and test. Each defaults to a single
#'   value matching [run_trial()]'s default.
#' @param n_chr,seg_sites,n_boot Held fixed across the grid; passed to
#'   [run_trial()].
#' @param n_rep Number of independent replicates per grid cell. Default `10`.
#' @param seed Optional integer seed; each replicate gets `seed + rep_index`
#'   for reproducibility without repeating the same simulation.
#'
#' @return A data frame stacking one [run_trial()] row per replicate, with
#'   an added `rep` column.
#' @export
benchmark_grid <- function(divergence = 1000, error_rate = 0, missing_rate = 0,
                            n_per_species = 5, n_chr = 10, seg_sites = 200,
                            n_boot = 500, n_rep = 10, seed = NULL) {
  grid <- expand.grid(
    divergence = divergence, error_rate = error_rate,
    missing_rate = missing_rate, n_per_species = n_per_species,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  rows <- list()
  for (i in seq_len(nrow(grid))) {
    for (r in seq_len(n_rep)) {
      trial_seed <- if (is.null(seed)) NULL else seed + (i - 1) * n_rep + r
      res <- run_trial(
        n_per_species = grid$n_per_species[i], divergence = grid$divergence[i],
        n_chr = n_chr, seg_sites = seg_sites,
        error_rate = grid$error_rate[i], missing_rate = grid$missing_rate[i],
        n_boot = n_boot, seed = trial_seed
      )
      res$rep <- r
      rows[[length(rows) + 1]] <- res
    }
  }

  do.call(rbind, rows)
}
