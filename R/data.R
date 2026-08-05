#' Construct an xpar genotype dataset
#'
#' Wraps a dosage genotype matrix together with the per-SNP linkage-group
#' (or scaffold/chromosome) labels that the block bootstrap resamples over.
#' Resampling by linkage group rather than by SNP is required because SNPs
#' within a group are not independent; treating them as independent would
#' understate uncertainty.
#'
#' @param geno A numeric matrix of allele dosages, samples in rows, SNPs in
#'   columns. Values must be `0`, `1`, `2`, or `NA`. Row names must be sample
#'   IDs; column names must be SNP IDs.
#' @param scaffold A character or factor vector, one entry per SNP (i.e.
#'   `length(scaffold) == ncol(geno)`), giving the linkage group, scaffold,
#'   or chromosome each SNP belongs to.
#'
#' @return An object of class `xpar_data`.
#' @export
#'
#' @examples
#' geno <- matrix(sample(c(0, 1, 2, NA), 40, replace = TRUE), nrow = 4,
#'                 dimnames = list(paste0("ind", 1:4), paste0("snp", 1:10)))
#' scaffold <- rep(paste0("chr", 1:2), each = 5)
#' xpar_data(geno, scaffold)
xpar_data <- function(geno, scaffold) {
  if (!is.matrix(geno) || !is.numeric(geno)) {
    stop("`geno` must be a numeric matrix.", call. = FALSE)
  }
  if (is.null(rownames(geno)) || is.null(colnames(geno))) {
    stop("`geno` must have row names (sample IDs) and column names (SNP IDs).",
         call. = FALSE)
  }
  bad_vals <- setdiff(unique(as.vector(geno)), c(0, 1, 2, NA))
  if (length(bad_vals) > 0) {
    stop("`geno` must contain only dosage values 0, 1, 2, or NA. Found: ",
         paste(utils::head(bad_vals, 5), collapse = ", "), call. = FALSE)
  }
  if (length(scaffold) != ncol(geno)) {
    stop("`scaffold` must have one entry per SNP (ncol(geno) = ", ncol(geno),
         ", length(scaffold) = ", length(scaffold), ").", call. = FALSE)
  }

  structure(
    list(
      geno = geno,
      scaffold = as.character(scaffold),
      sample_id = rownames(geno),
      snp_id = colnames(geno)
    ),
    class = "xpar_data"
  )
}

#' @export
print.xpar_data <- function(x, ...) {
  cat(sprintf(
    "<xpar_data> %d samples x %d SNPs across %d linkage groups\n",
    length(x$sample_id), length(x$snp_id), length(unique(x$scaffold))
  ))
  invisible(x)
}

check_ids <- function(data, ids, arg_name) {
  missing <- setdiff(ids, data$sample_id)
  if (length(missing) > 0) {
    stop("`", arg_name, "` not found in data: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}
