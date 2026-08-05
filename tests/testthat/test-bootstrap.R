test_that("block_bootstrap identifies the true pair with win_prob 1", {
  d <- make_scan_fixture()
  boot <- block_bootstrap(d, candidates = c("P1", "P2", "P3"), offspring = "Off",
                           self = FALSE, n_boot = 300, seed = 1)

  expect_s3_class(boot, "xpar_boot")
  top <- rank_pairs(boot, offspring = "Off", top_n = 1)
  expect_equal(sort(c(top$parent1, top$parent2)), c("P1", "P2"))
  expect_equal(top$rate, 0)
  # P1xP2 has zero errors in every linkage group, so it is tied for the
  # bootstrap minimum in every replicate and (being first in pair order)
  # always wins the tie-break -> win probability is exactly 1.
  expect_equal(top$win_prob, 1)
})

test_that("block_bootstrap is reproducible with the same seed", {
  d <- make_scan_fixture()
  b1 <- block_bootstrap(d, c("P1", "P2", "P3"), "Off", self = FALSE, n_boot = 100, seed = 42)
  b2 <- block_bootstrap(d, c("P1", "P2", "P3"), "Off", self = FALSE, n_boot = 100, seed = 42)
  expect_equal(b1$summary$rate, b2$summary$rate)
  expect_equal(b1$summary$win_prob, b2$summary$win_prob)
})

test_that("block_bootstrap requires at least 2 pairs", {
  d <- make_scan_fixture()
  expect_error(
    block_bootstrap(d, candidates = "P1", offspring = "Off", self = TRUE, n_boot = 10),
    "at least 2"
  )
})

test_that("rank_pairs errors on unknown offspring", {
  d <- make_scan_fixture()
  boot <- block_bootstrap(d, c("P1", "P2", "P3"), "Off", self = FALSE, n_boot = 20, seed = 1)
  expect_error(rank_pairs(boot, offspring = "NotThere"), "not found")
})
