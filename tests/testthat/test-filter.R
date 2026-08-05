test_that("filter_by_group restricts to allowed group combinations and re-ranks", {
  d <- make_scan_fixture()
  boot <- block_bootstrap(d, c("P1", "P2", "P3"), "Off", self = FALSE, n_boot = 200, seed = 7)

  group <- c(P1 = "sp_a", P2 = "sp_b", P3 = "sp_a")
  allowed <- data.frame(g1 = "sp_a", g2 = "sp_b", stringsAsFactors = FALSE)

  filtered <- filter_by_group(boot, group, allowed)
  s <- filtered$by_offspring[["Off"]]$summary

  # Only P1xP2 (sp_a x sp_b) and P2xP3 (sp_b x sp_a) are allowed;
  # P1xP3 (sp_a x sp_a) must be dropped.
  expect_equal(nrow(s), 2)
  pair_sets <- lapply(seq_len(nrow(s)), function(i) sort(c(s$parent1[i], s$parent2[i])))
  expect_false(any(vapply(pair_sets, function(p) identical(p, c("P1", "P3")), logical(1))))
  expect_equal(sum(s$win_prob), 1)
})

test_that("filter_by_group with allowed = NULL keeps all pairs but adds labels", {
  d <- make_scan_fixture()
  boot <- block_bootstrap(d, c("P1", "P2", "P3"), "Off", self = FALSE, n_boot = 50, seed = 3)
  group <- c(P1 = "sp_a", P2 = "sp_b", P3 = "sp_a")

  out <- filter_by_group(boot, group)
  s <- out$by_offspring[["Off"]]$summary
  expect_equal(nrow(s), 3)
  expect_true(all(c("group1", "group2") %in% names(s)))
})
