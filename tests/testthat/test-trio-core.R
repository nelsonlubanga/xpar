test_that("trio_consistency implements the full Mendelian rule table", {
  # One valid and one invalid offspring genotype for each of the 6 unordered
  # parent-genotype-class combinations, plus one locus missing in the
  # offspring and one missing in a parent.
  p1  <- c(0, 0,   0, 0,   0, 0,   1, 1,   1, 1,   2, 2,   0,  NA)
  p2  <- c(0, 0,   1, 1,   2, 2,   1, 1,   2, 2,   2, 2,   0,  0)
  off <- c(0, 1,   1, 2,   1, 0,   2, 0,   2, 0,   2, 1,   NA, 0)
  expect_error <- c(FALSE, TRUE,  FALSE, TRUE,  FALSE, TRUE,
                     FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE,
                     FALSE, FALSE)
  expect_called <- c(rep(TRUE, 12), FALSE, FALSE)

  geno <- rbind(P1 = p1, P2 = p2, Off = off)
  colnames(geno) <- paste0("snp", seq_along(p1))
  storage.mode(geno) <- "numeric"
  d <- xpar_data(geno, rep("chr1", length(p1)))

  res <- trio_consistency(d, "P1", "P2", "Off")
  expect_equal(res$called, expect_called, ignore_attr = "names")
  expect_equal(res$error, expect_error & expect_called, ignore_attr = "names")

  er <- error_rate(res)
  expect_equal(er$n_called, 12)
  expect_equal(er$n_error, 5)
  expect_equal(er$rate, 5 / 12)
})

test_that("trio_consistency rejects unknown sample IDs", {
  geno <- matrix(0, nrow = 2, ncol = 2, dimnames = list(c("A", "B"), c("s1", "s2")))
  d <- xpar_data(geno, c("chr1", "chr1"))
  expect_error(trio_consistency(d, "A", "B", "nope"), "not found")
})

test_that("mendelian_error flags only opposing-homozygote loci", {
  parent <- c(0, 0, 1, 2, 2, NA)
  off    <- c(0, 2, 1, 0, 2, 1)
  geno <- rbind(Pa = parent, Off = off)
  colnames(geno) <- paste0("snp", 1:6)
  d <- xpar_data(geno, rep("chr1", 6))

  res <- mendelian_error(d, "Pa", "Off")
  expect_equal(res$called, c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE), ignore_attr = "names")
  expect_equal(res$error,  c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE), ignore_attr = "names")

  er <- error_rate(res)
  expect_equal(er$n_called, 5)
  expect_equal(er$n_error, 2)
  expect_equal(er$rate, 0.4)
})

test_that("trio_scan ranks the true pair above decoys", {
  d <- make_scan_fixture()
  out <- trio_scan(d, candidates = c("P1", "P2", "P3"), offspring = "Off", self = FALSE)

  expect_equal(nrow(out), 3)
  best <- out[1, ]
  expect_equal(sort(c(best$parent1, best$parent2)), c("P1", "P2"))
  expect_equal(best$rate, 0)
})
