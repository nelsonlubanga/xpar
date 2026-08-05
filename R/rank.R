#' Rank candidate pairs by incompatibility rate
#'
#' Returns a rate-sorted view of the results in an `xpar_boot` object.
#'
#' @param x An `xpar_boot` object, as returned by [block_bootstrap()].
#' @param offspring Optional character vector restricting the output to
#'   specific offspring. Default (`NULL`) returns all offspring.
#' @param top_n Optional integer; keep only the top `top_n` pairs (by
#'   ascending rate) per offspring.
#'
#' @return A data frame sorted by `offspring`, then ascending `rate`.
#' @export
rank_pairs <- function(x, offspring = NULL, top_n = NULL) {
  stopifnot(inherits(x, "xpar_boot"))
  offs <- if (is.null(offspring)) names(x$by_offspring) else offspring
  missing_off <- setdiff(offs, names(x$by_offspring))
  if (length(missing_off) > 0) {
    stop("Offspring not found in `x`: ", paste(missing_off, collapse = ", "),
         call. = FALSE)
  }

  out <- do.call(rbind, lapply(offs, function(off) {
    s <- x$by_offspring[[off]]$summary
    s <- s[order(s$rate), ]
    if (!is.null(top_n)) s <- utils::head(s, top_n)
    s
  }))
  rownames(out) <- NULL
  out
}
