#' Per-locus Mendelian trio consistency
#'
#' For each SNP, determines whether the offspring's genotype is a possible
#' product of the two candidate parents' genotypes under diploid Mendelian
#' inheritance (i.e. bi-allelic dosage 0/1/2 transmission), and flags
#' incompatible (impossible) combinations.
#'
#' @param data An `xpar_data` object.
#' @param parent1,parent2,offspring Sample IDs present in `data`.
#'
#' @return A list with logical vectors `called` (locus had non-missing
#'   genotypes for all three individuals) and `error` (locus was called and
#'   the offspring genotype was impossible given the two parents), both of
#'   length `ncol(data$geno)`.
#' @export
trio_consistency <- function(data, parent1, parent2, offspring) {
  check_ids(data, c(parent1, parent2, offspring), "parent1/parent2/offspring")

  p1 <- data$geno[parent1, ]
  p2 <- data$geno[parent2, ]
  off <- data$geno[offspring, ]

  called <- !is.na(p1) & !is.na(p2) & !is.na(off)
  possible <- rep(FALSE, length(off))

  is_hom_ref <- called & p1 == 0 & p2 == 0
  possible[is_hom_ref] <- off[is_hom_ref] == 0

  is_ref_het <- called & ((p1 == 0 & p2 == 1) | (p1 == 1 & p2 == 0))
  possible[is_ref_het] <- off[is_ref_het] %in% c(0, 1)

  is_ref_alt <- called & ((p1 == 0 & p2 == 2) | (p1 == 2 & p2 == 0))
  possible[is_ref_alt] <- off[is_ref_alt] == 1

  is_het_het <- called & p1 == 1 & p2 == 1
  possible[is_het_het] <- TRUE

  is_het_alt <- called & ((p1 == 1 & p2 == 2) | (p1 == 2 & p2 == 1))
  possible[is_het_alt] <- off[is_het_alt] %in% c(1, 2)

  is_hom_alt <- called & p1 == 2 & p2 == 2
  possible[is_hom_alt] <- off[is_hom_alt] == 2

  list(called = called, error = called & !possible)
}

#' Single-candidate parent-offspring incompatibility
#'
#' Screens one candidate parent against one offspring using the
#' opposing-homozygote rule: a locus is incompatible if parent and offspring
#' are both homozygous for opposite alleles (0 vs 2), which is impossible
#' under Mendelian transmission regardless of the other parent's genotype.
#' This is a fast triage step, cheaper but less specific than
#' [trio_consistency()].
#'
#' @inheritParams trio_consistency
#' @param parent,offspring Sample IDs present in `data`.
#'
#' @return A list with logical vectors `called` and `error`.
#' @export
mendelian_error <- function(data, parent, offspring) {
  check_ids(data, c(parent, offspring), "parent/offspring")

  p <- data$geno[parent, ]
  off <- data$geno[offspring, ]

  called <- !is.na(p) & !is.na(off)
  error <- called & ((p == 0 & off == 2) | (p == 2 & off == 0))

  list(called = called, error = error)
}

#' Summarize a called/error locus list into a rate
#'
#' @param x A list as returned by [trio_consistency()] or [mendelian_error()].
#' @return A one-row data frame with `n_called`, `n_error`, and `rate`.
#' @export
error_rate <- function(x) {
  n_called <- sum(x$called)
  n_error <- sum(x$error)
  data.frame(
    n_called = n_called,
    n_error = n_error,
    rate = if (n_called > 0) n_error / n_called else NA_real_
  )
}

#' Scan a candidate pool for the best-fitting parent pair(s)
#'
#' Runs [trio_consistency()] over every candidate pair (or a supplied subset
#' of pairs) against one or more offspring, and returns a tidy ranking by
#' incompatibility rate. This is the core "blind search" driver: no parent
#' needs to be assumed known in advance.
#'
#' @param data An `xpar_data` object.
#' @param candidates Character vector of candidate parent sample IDs.
#' @param offspring Character vector of offspring sample IDs to test.
#' @param pairs Optional two-column character matrix or data frame of
#'   specific `(parent1, parent2)` pairs to restrict the scan to. If `NULL`
#'   (the default), all unique unordered pairs of `candidates` are tested.
#' @param self Logical; if `TRUE` (default), also test each candidate against
#'   itself (a self-cross / selfing pair).
#'
#' @return A data frame with one row per (offspring, parent1, parent2) triad,
#'   including `n_called`, `n_error`, and `rate`, sorted by `offspring` then
#'   ascending `rate`.
#' @export
trio_scan <- function(data, candidates, offspring, pairs = NULL, self = TRUE) {
  check_ids(data, candidates, "candidates")
  check_ids(data, offspring, "offspring")

  if (is.null(pairs)) {
    pairs <- if (length(candidates) >= 2) utils::combn(candidates, 2, simplify = FALSE) else list()
    if (self) {
      pairs <- c(pairs, lapply(candidates, function(x) c(x, x)))
    }
    if (length(pairs) == 0) {
      stop("No candidate pairs to test: provide >= 2 candidates, or self = TRUE.",
           call. = FALSE)
    }
  } else {
    pairs <- lapply(seq_len(nrow(pairs)), function(i) as.character(pairs[i, 1:2]))
  }

  rows <- lapply(offspring, function(off) {
    do.call(rbind, lapply(pairs, function(pr) {
      res <- trio_consistency(data, pr[1], pr[2], off)
      cbind(
        data.frame(offspring = off, parent1 = pr[1], parent2 = pr[2],
                    stringsAsFactors = FALSE),
        error_rate(res)
      )
    }))
  })

  out <- do.call(rbind, rows)
  out[order(out$offspring, out$rate), ]
}
