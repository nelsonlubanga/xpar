test_that("calibrate_reference measures the rate for a trusted trio", {
  d <- make_scan_fixture()
  ref <- calibrate_reference(d, "P1", "P2", "Off")

  expect_equal(ref$parent1, "P1")
  expect_equal(ref$parent2, "P2")
  expect_equal(ref$offspring, "Off")
  expect_equal(ref$n_called, 12)
  expect_equal(ref$n_error, 0)
  expect_equal(ref$rate, 0)
})

test_that("test_vs_reference flags a candidate that departs from the reference", {
  d <- make_scan_fixture()
  ref <- calibrate_reference(d, "P1", "P2", "Off")
  candidate <- error_rate(trio_consistency(d, "P1", "P3", "Off"))
  candidate$parent1 <- "P1"; candidate$parent2 <- "P3"

  # Small counts trigger prop.test's routine chi-squared-approximation
  # warning; expected here, not a bug.
  ht <- suppressWarnings(test_vs_reference(candidate, ref))
  expect_s3_class(ht, "htest")
  expect_true(is.numeric(ht$p.value))
})
