#' Re-rank candidate pairs using independent group/species labels
#'
#' Restricts the candidate pairs in an `xpar_boot` object to those
#' consistent with independently established group membership (e.g. species,
#' subspecies, or known population of origin), then recomputes win
#' probabilities *within that restricted subset* using the original
#' bootstrap replicates. This is a post-hoc consistency check, not a search
#' seed: [block_bootstrap()] should always be run first on the full,
#' unconstrained candidate pool so that the blind result can be compared
#' against the group-constrained one.
#'
#' @param x An `xpar_boot` object, as returned by [block_bootstrap()].
#' @param group A named character vector mapping sample ID to group label
#'   (e.g. `c(candidate1 = "sp_a", candidate2 = "sp_b")`).
#' @param allowed Optional two-column character matrix or data frame of
#'   allowed `(group1, group2)` combinations (order-insensitive). If `NULL`
#'   (the default), every pair present in `x` is kept and only `group1` /
#'   `group2` labels are attached — use this to inspect group composition
#'   before deciding what to allow.
#'
#' @return An `xpar_boot` object with `group1`/`group2` columns added to each
#'   summary, restricted to `allowed` combinations (if supplied), and
#'   `win_prob` recomputed within the restricted subset.
#' @export
filter_by_group <- function(x, group, allowed = NULL) {
  stopifnot(inherits(x, "xpar_boot"))

  group_of <- function(id) unname(group[id])
  pair_allowed <- function(a, b) {
    if (is.null(allowed)) return(TRUE)
    any((allowed[, 1] == a & allowed[, 2] == b) |
        (allowed[, 1] == b & allowed[, 2] == a))
  }

  for (off in names(x$by_offspring)) {
    entry <- x$by_offspring[[off]]
    s <- entry$summary
    s$group1 <- group_of(s$parent1)
    s$group2 <- group_of(s$parent2)
    keep <- mapply(pair_allowed, s$group1, s$group2)

    entry$summary <- s[keep, , drop = FALSE]
    entry$boot_rate <- entry$boot_rate[keep, , drop = FALSE]
    entry$pairs <- entry$pairs[keep]

    if (nrow(entry$summary) > 0) {
      win <- apply(entry$boot_rate, 2, function(col) {
        if (all(is.na(col))) return(NA_integer_)
        which.min(col)
      })
      entry$summary$win_prob <- tabulate(win, nbins = nrow(entry$summary)) / x$n_boot
    }

    x$by_offspring[[off]] <- entry
  }

  x$summary <- do.call(rbind, lapply(x$by_offspring, `[[`, "summary"))
  x
}
