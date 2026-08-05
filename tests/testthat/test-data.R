test_that("xpar_data validates inputs", {
  geno <- matrix(c(0, 1, 2, NA), nrow = 1, dimnames = list("s1", paste0("snp", 1:4)))

  expect_s3_class(xpar_data(geno, c("c1", "c1", "c2", "c2")), "xpar_data")

  expect_error(xpar_data(as.data.frame(geno), c("c1", "c1", "c2", "c2")),
               "numeric matrix")

  geno_no_names <- geno
  dimnames(geno_no_names) <- NULL
  expect_error(xpar_data(geno_no_names, c("c1", "c1", "c2", "c2")),
               "row names")

  bad_geno <- geno
  bad_geno[1, 1] <- 3
  expect_error(xpar_data(bad_geno, c("c1", "c1", "c2", "c2")),
               "dosage values")

  expect_error(xpar_data(geno, c("c1", "c2")), "one entry per SNP")
})

test_that("print.xpar_data reports dimensions", {
  geno <- matrix(0, nrow = 2, ncol = 4,
                  dimnames = list(c("s1", "s2"), paste0("snp", 1:4)))
  d <- xpar_data(geno, rep(c("chr1", "chr2"), each = 2))
  expect_output(print(d), "2 samples x 4 SNPs across 2 linkage groups")
})
