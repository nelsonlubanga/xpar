#' Block bootstrap over linkage groups
#'
#' Quantifies uncertainty in candidate-pair incompatibility rates by
#' resampling whole linkage groups (scaffolds/chromosomes) with replacement,
#' not individual SNPs. SNPs within a linkage group are correlated by
#' linkage disequilibrium, so resampling at the SNP level would understate
#' uncertainty; resampling at the linkage-group level treats each group as
#' one (approximately independent) unit.
#'
#' For each bootstrap replicate, every candidate pair's incompatibility rate
#' is recomputed from the resampled linkage groups, and the pair with the
#' strict lowest rate in that replicate is recorded as the "winner". The
#' `win_prob` reported per pair is the fraction of replicates in which it
#' won — a direct, assumption-light measure of how often a pair is the best
#' explanation for the offspring's genotype relative to every other
#' candidate pair tested, given the uncertainty in the data.
#'
#' @inheritParams trio_scan
#' @param n_boot Number of bootstrap replicates. Default `1000`.
#' @param conf_level Confidence level for the percentile interval on each
#'   pair's rate. Default `0.95`.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `xpar_boot`: a list with `summary` (all
#'   offspring stacked), `by_offspring` (per-offspring summary plus the raw
#'   `n_pairs x n_boot` bootstrap rate matrix), `n_boot`, and `conf_level`.
#' @export
block_bootstrap <- function(data, candidates, offspring, pairs = NULL,
                             self = TRUE, n_boot = 1000, conf_level = 0.95,
                             seed = NULL) {
  check_ids(data, candidates, "candidates")
  check_ids(data, offspring, "offspring")
  if (!is.null(seed)) set.seed(seed)

  if (is.null(pairs)) {
    n_unordered <- choose(length(candidates), 2) + if (self) length(candidates) else 0
    if (n_unordered < 2) {
      stop("Need at least 2 candidate pairs to compute win probabilities ",
           "(got ", length(candidates), " candidate(s), self = ", self, ").",
           call. = FALSE)
    }
    pairs <- if (length(candidates) >= 2) utils::combn(candidates, 2, simplify = FALSE) else list()
    if (self) pairs <- c(pairs, lapply(candidates, function(x) c(x, x)))
  } else {
    pairs <- lapply(seq_len(nrow(pairs)), function(i) as.character(pairs[i, 1:2]))
  }
  n_pairs <- length(pairs)
  if (n_pairs < 2) {
    stop("Need at least 2 candidate pairs to compute win probabilities.",
         call. = FALSE)
  }

  scaffold_f <- factor(data$scaffold, levels = unique(data$scaffold))
  n_scaffold <- nlevels(scaffold_f)
  alpha <- 1 - conf_level

  results <- vector("list", length(offspring))
  names(results) <- offspring

  for (off in offspring) {
    called_mat <- matrix(0, n_scaffold, n_pairs)
    error_mat  <- matrix(0, n_scaffold, n_pairs)
    for (j in seq_len(n_pairs)) {
      res <- trio_consistency(data, pairs[[j]][1], pairs[[j]][2], off)
      called_mat[, j] <- as.numeric(rowsum(as.numeric(res$called), scaffold_f, reorder = TRUE))
      error_mat[, j]  <- as.numeric(rowsum(as.numeric(res$error),  scaffold_f, reorder = TRUE))
    }

    boot_idx <- replicate(n_boot, tabulate(
      sample.int(n_scaffold, n_scaffold, replace = TRUE), nbins = n_scaffold
    ))
    boot_called <- t(called_mat) %*% boot_idx  # n_pairs x n_boot
    boot_error  <- t(error_mat)  %*% boot_idx
    boot_rate   <- boot_error / boot_called

    win <- apply(boot_rate, 2, function(col) {
      if (all(is.na(col))) return(NA_integer_)
      which.min(col)
    })
    win_prob <- tabulate(win, nbins = n_pairs) / n_boot

    point_called <- colSums(called_mat)
    point_error  <- colSums(error_mat)
    point_rate   <- point_error / point_called
    ci <- t(apply(boot_rate, 1, stats::quantile,
                  probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE))

    summary_df <- data.frame(
      offspring = off,
      parent1 = vapply(pairs, `[`, character(1), 1),
      parent2 = vapply(pairs, `[`, character(1), 2),
      n_called = point_called,
      n_error = point_error,
      rate = point_rate,
      ci_low = ci[, 1],
      ci_high = ci[, 2],
      win_prob = win_prob,
      stringsAsFactors = FALSE
    )
    # Kept in `pairs` order (not sorted by rate) so that row i of `summary`,
    # row i of `boot_rate`, and pairs[[i]] always refer to the same pair —
    # this is what lets filter_by_group() subset all three consistently.
    # Use rank_pairs() to get a rate-sorted view.

    results[[off]] <- list(summary = summary_df, boot_rate = boot_rate, pairs = pairs)
  }

  structure(
    list(
      by_offspring = results,
      summary = do.call(rbind, lapply(results, `[[`, "summary")),
      n_boot = n_boot,
      conf_level = conf_level
    ),
    class = "xpar_boot"
  )
}

#' @export
print.xpar_boot <- function(x, ...) {
  cat(sprintf("<xpar_boot> %d bootstrap replicates, %d%% CI\n",
              x$n_boot, round(x$conf_level * 100)))
  for (off in names(x$by_offspring)) {
    cat(sprintf("\n-- %s --\n", off))
    top <- utils::head(x$by_offspring[[off]]$summary, 3)
    print(top, row.names = FALSE)
  }
  invisible(x)
}
