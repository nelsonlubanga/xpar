test_that("simulate_hybrid_pedigree produces a consistent, recoverable trio", {
  testthat::skip_if_not_installed("AlphaSimR")

  sim <- simulate_hybrid_pedigree(n_per_species = 3, divergence = 2000,
                                   n_chr = 4, seg_sites = 100, seed = 1)

  expect_equal(dim(sim$geno)[1], 7)  # 3 + 3 candidates + 1 offspring
  expect_true(all(c(sim$true_parent1, sim$true_parent2, sim$offspring) %in% rownames(sim$geno)))
  expect_equal(length(sim$candidates), 6)
  expect_equal(length(sim$scaffold), ncol(sim$geno))
  expect_true(all(sim$species[sim$candidates] %in% c("spA", "spB")))

  # The true parents should be genuinely Mendelian-consistent with the
  # offspring: at low divergence / zero injected noise, error rate should
  # be low (real recombination/mutation noise aside, not literally zero).
  d <- xpar_data(sim$geno, sim$scaffold)
  res <- trio_consistency(d, sim$true_parent1, sim$true_parent2, sim$offspring)
  er <- error_rate(res)
  expect_lt(er$rate, 0.05)
})

test_that("simulate_hybrid_pedigree errors clearly without AlphaSimR", {
  testthat::skip_if(requireNamespace("AlphaSimR", quietly = TRUE),
                     "AlphaSimR is installed; can't test the missing-dependency path")
  expect_error(simulate_hybrid_pedigree(), "AlphaSimR")
})

test_that("inject_noise adds the requested amount of error and missingness", {
  set.seed(1)
  geno <- matrix(sample(0:2, 20000, replace = TRUE), nrow = 20)

  noisy <- inject_noise(geno, error_rate = 0.1, missing_rate = 0.2, seed = 2)
  frac_missing <- mean(is.na(noisy))
  expect_equal(frac_missing, 0.2, tolerance = 0.02)

  frac_changed <- mean(noisy[!is.na(noisy)] != geno[!is.na(noisy)])
  # error and missingness are drawn independently, so conditioning on "not
  # missing" doesn't change the error-flip probability.
  expect_equal(frac_changed, 0.1, tolerance = 0.02)

  expect_true(all(noisy %in% c(0, 1, 2, NA)))
})

test_that("inject_noise with all-zero rates returns the input unchanged", {
  geno <- matrix(c(0, 1, 2, NA, 1, 0), nrow = 2)
  expect_equal(inject_noise(geno, 0, 0), geno)
})

test_that("run_trial returns a well-formed one-row result", {
  testthat::skip_if_not_installed("AlphaSimR")

  res <- run_trial(n_per_species = 3, divergence = 2000, n_chr = 4, seg_sites = 100,
                    error_rate = 0, missing_rate = 0, n_boot = 50, seed = 1)

  expect_equal(nrow(res), 1)
  expect_true(res$true_pair_rank >= 1)
  expect_true(is.logical(res$true_pair_is_top))
  expect_true(res$true_pair_win_prob >= 0 && res$true_pair_win_prob <= 1)
})

test_that("benchmark_grid stacks one row per condition x replicate", {
  testthat::skip_if_not_installed("AlphaSimR")

  out <- benchmark_grid(
    divergence = c(500, 2000), error_rate = 0, missing_rate = 0,
    n_per_species = 3, n_chr = 4, seg_sites = 100, n_boot = 50,
    n_rep = 2, seed = 10
  )

  expect_equal(nrow(out), 2 * 2)  # 2 divergence values x 2 reps
  expect_true(all(c("divergence", "rep", "true_pair_is_top") %in% names(out)))
})
