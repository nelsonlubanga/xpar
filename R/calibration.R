#' Empirical reference (self-calibration) rate
#'
#' The core departure from tools that assume a fixed genotyping-error rate
#' (e.g. sequoia's `Err` parameter): instead of assuming what the
#' "background" incompatibility rate should be, `calibrate_reference()`
#' measures it directly from a trio you already trust — either a pedigree-
#' confirmed true parent pair, or (in a fully blind analysis) the
#' best-supported pair from [trio_scan()], clearly labeled as such. Every
#' other candidate pair is then judged against this data-derived reference
#' rather than a universal constant, which is what lets the same method work
#' whether the true discordance is dominated by genotyping error (a single
#' population) or by real allele-frequency divergence (an interspecific
#' cross).
#'
#' @inheritParams trio_consistency
#'
#' @return A one-row data frame: `parent1`, `parent2`, `offspring`,
#'   `n_called`, `n_error`, `rate`.
#' @export
calibrate_reference <- function(data, parent1, parent2, offspring) {
  res <- trio_consistency(data, parent1, parent2, offspring)
  cbind(
    data.frame(parent1 = parent1, parent2 = parent2, offspring = offspring,
                stringsAsFactors = FALSE),
    error_rate(res)
  )
}

#' Test a candidate pair's rate against the empirical reference
#'
#' Two-proportion test (via [stats::prop.test()]) comparing a candidate
#' pair's incompatibility rate to the empirical reference rate produced by
#' [calibrate_reference()]. A candidate pair statistically indistinguishable
#' from the reference is consistent with being a true parent pair; a rate
#' significantly higher than the reference is not.
#'
#' @param candidate A one-row data frame as returned by a single row of
#'   [trio_scan()] (must contain `n_error` and `n_called`).
#' @param reference A one-row data frame as returned by
#'   [calibrate_reference()] (must contain `n_error` and `n_called`).
#' @param ... Passed to [stats::prop.test()].
#'
#' @return The `htest` object returned by [stats::prop.test()].
#' @export
test_vs_reference <- function(candidate, reference, ...) {
  stats::prop.test(
    x = c(candidate$n_error, reference$n_error),
    n = c(candidate$n_called, reference$n_called),
    ...
  )
}
