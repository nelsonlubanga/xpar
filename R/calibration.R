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
#' When both the candidate and the reference have zero errors,
#' [stats::prop.test()] degenerates (`NaN` statistic, `NA` p-value, with a
#' spurious "Chi-squared approximation may be incorrect" warning) because
#' the pooled proportion is exactly 0. That case — a candidate matching the
#' reference perfectly — is unambiguous evidence of no difference, so it is
#' short-circuited to a p-value of `1` instead.
#'
#' @param candidate A one-row data frame as returned by a single row of
#'   [trio_scan()] (must contain `n_error` and `n_called`).
#' @param reference A one-row data frame as returned by
#'   [calibrate_reference()] (must contain `n_error` and `n_called`).
#' @param ... Passed to [stats::prop.test()]. Ignored in the zero-vs-zero
#'   case, except for `conf.level` (default `0.95`).
#'
#' @return An `htest` object, as returned by [stats::prop.test()] (or an
#'   equivalent stub in the zero-vs-zero case).
#' @export
test_vs_reference <- function(candidate, reference, ...) {
  if (candidate$n_error == 0 && reference$n_error == 0) {
    return(zero_vs_zero_prop_test(...))
  }
  stats::prop.test(
    x = c(candidate$n_error, reference$n_error),
    n = c(candidate$n_called, reference$n_called),
    ...
  )
}

#' Stub `htest` for a candidate and reference that both have zero errors
#'
#' @param ... Only `conf.level` is used (default `0.95`); other arguments
#'   (e.g. `alternative`, `correct`) are accepted and ignored for
#'   compatibility with [stats::prop.test()] call sites.
#' @return An `htest` object with `p.value = 1`.
#' @noRd
zero_vs_zero_prop_test <- function(...) {
  dots <- list(...)
  conf_level <- if (is.null(dots$conf.level)) 0.95 else dots$conf.level

  structure(
    list(
      statistic = c("X-squared" = 0),
      parameter = c(df = 1),
      p.value = 1,
      estimate = c(`prop 1` = 0, `prop 2` = 0),
      null.value = c(`difference in proportions` = 0),
      conf.int = structure(c(0, 0), conf.level = conf_level),
      alternative = "two.sided",
      method = paste(
        "2-sample test for equality of proportions",
        "(degenerate: both rates are exactly 0)"
      ),
      data.name = "candidate and reference"
    ),
    class = "htest"
  )
}
