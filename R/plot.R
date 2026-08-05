#' Plot ranked candidate pairs
#'
#' Bar/point plot of candidate-pair incompatibility rates with bootstrap
#' confidence intervals, faceted by offspring. Requires the \pkg{ggplot2}
#' package.
#'
#' @param x An `xpar_boot` object, as returned by [block_bootstrap()].
#' @param top_n Number of top (lowest-rate) pairs to show per offspring.
#'   Default `10`.
#' @param highlight Optional character vector of length 2 giving a
#'   `(parent1, parent2)` pair to outline, e.g. a pedigree-documented pair,
#'   for visual comparison against the data-driven ranking.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#' @export
plot.xpar_boot <- function(x, top_n = 10, highlight = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" is required for plot.xpar_boot(). ",
         "Install it with install.packages(\"ggplot2\").", call. = FALSE)
  }

  df <- rank_pairs(x, top_n = top_n)
  df$pair_label <- paste(df$parent1, "+", df$parent2)
  df$pair_label <- stats::reorder(df$pair_label, -df$rate)

  if (is.null(highlight)) {
    df$highlighted <- FALSE
  } else {
    df$highlighted <- (df$parent1 == highlight[1] & df$parent2 == highlight[2]) |
      (df$parent1 == highlight[2] & df$parent2 == highlight[1])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = rate, y = pair_label)) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      height = 0, color = "#2a78d6"
    ) +
    ggplot2::geom_point(size = 3, color = "#2a78d6") +
    ggplot2::facet_wrap(~offspring, scales = "free_y") +
    ggplot2::labs(
      x = "Mendelian incompatibility rate (lower = more parent-like)",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(highlight) && any(df$highlighted)) {
    p <- p + ggplot2::geom_point(
      data = df[df$highlighted, , drop = FALSE],
      shape = 21, size = 6.5, stroke = 1.1, color = "#0b0b0b", fill = NA
    )
  }

  p
}
