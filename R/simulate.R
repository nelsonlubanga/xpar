#' Simulate a two-species candidate pool and a known interspecific cross
#'
#' Wraps \pkg{AlphaSimR}'s coalescent founder simulation (`runMacs()`) to
#' generate two diverged founder pools — using the `split` argument to set
#' how many generations ago they shared a common ancestor — then crosses one
#' individual from each to produce a known F1 hybrid offspring. The
#' remaining individuals in both pools serve as decoy candidates, giving a
#' candidate-pool parentage problem with a fully known ground truth, at a
#' controllable level of interspecific divergence.
#'
#' Requires the \pkg{AlphaSimR} package.
#'
#' @param n_per_species Number of individuals simulated per species pool
#'   (one from each side is used as the true parent; the rest are decoys).
#'   Default `5`.
#' @param divergence Generations since the two species' common ancestor,
#'   passed to `AlphaSimR::runMacs(split = divergence)`. Larger values mean
#'   more diverged species (more locus-heterogeneous "error" relative to a
#'   single-population model). Default `1000`.
#' @param n_chr Number of chromosomes to simulate. Default `10`.
#' @param seg_sites Segregating sites retained per chromosome. Default `200`.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with:
#'   \item{geno}{Dosage matrix (0/1/2), samples x SNPs, ready for
#'     [xpar_data()].}
#'   \item{scaffold}{Per-SNP chromosome label, for the `scaffold` argument
#'     of [xpar_data()].}
#'   \item{candidates}{Character vector of all candidate parent IDs (true
#'     parents + decoys, both species).}
#'   \item{species}{Named character vector mapping candidate ID to species
#'     label (`"spA"` / `"spB"`), for [filter_by_group()].}
#'   \item{true_parent1, true_parent2}{IDs of the actual simulated parents.}
#'   \item{offspring}{ID of the simulated F1 hybrid.}
#' @export
simulate_hybrid_pedigree <- function(n_per_species = 5, divergence = 1000,
                                      n_chr = 10, seg_sites = 200, seed = NULL) {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("Package \"AlphaSimR\" is required for simulate_hybrid_pedigree(). ",
         "Install it with install.packages(\"AlphaSimR\").", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  n_total <- 2 * n_per_species
  founders <- AlphaSimR::runMacs(
    nInd = n_total, nChr = n_chr, segSites = seg_sites,
    species = "GENERIC", split = divergence
  )
  sp <- AlphaSimR::SimParam$new(founders)
  pop <- AlphaSimR::newPop(founders, simParam = sp)

  # runMacs(nInd, split = ...) returns the first half of individuals drawn
  # from one diverged lineage and the second half from the other.
  sp_a_id <- pop@id[seq_len(n_per_species)]
  sp_b_id <- pop@id[n_per_species + seq_len(n_per_species)]

  true_mother <- sp_a_id[1]
  true_father <- sp_b_id[1]
  cross_plan <- matrix(c(true_mother, true_father), nrow = 1)
  offspring_pop <- AlphaSimR::makeCross(pop, cross_plan, nProgeny = 1, simParam = sp)

  combined <- c(pop, offspring_pop)
  geno <- AlphaSimR::pullSegSiteGeno(combined, simParam = sp)

  id_map <- c(
    stats::setNames(paste0("spA_", seq_len(n_per_species)), sp_a_id),
    stats::setNames(paste0("spB_", seq_len(n_per_species)), sp_b_id),
    stats::setNames("Hybrid_off", offspring_pop@id)
  )
  rownames(geno) <- id_map[rownames(geno)]
  storage.mode(geno) <- "numeric"

  scaffold <- sub("_.*$", "", colnames(geno))

  candidates <- unname(id_map[c(sp_a_id, sp_b_id)])
  species <- stats::setNames(
    c(rep("spA", n_per_species), rep("spB", n_per_species)),
    candidates
  )

  list(
    geno = geno,
    scaffold = scaffold,
    candidates = candidates,
    species = species,
    true_parent1 = id_map[[true_mother]],
    true_parent2 = id_map[[true_father]],
    offspring = "Hybrid_off"
  )
}

#' Overlay genotyping error and missingness on a dosage matrix
#'
#' `AlphaSimR` returns true, error-free genotypes. Real data never is —
#' this overlays a controllable amount of genotyping error (a one-step
#' miscall: 0<->1 or 1<->2, reflecting a single misread allele, rather than
#' an arbitrary 0<->2 flip) and missingness, independently, so that
#' [simulate_hybrid_pedigree()] output can be used to test how a method's
#' power degrades as data quality worsens.
#'
#' @param geno A dosage matrix (0/1/2, optionally with `NA`).
#' @param error_rate Per-genotyped-call probability of a one-step miscall.
#'   Default `0`.
#' @param missing_rate Per-call probability of being set to `NA`, applied
#'   independently after `error_rate`. Default `0`.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A dosage matrix of the same dimensions as `geno`.
#' @export
inject_noise <- function(geno, error_rate = 0, missing_rate = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  out <- geno
  called <- which(!is.na(out))

  if (error_rate > 0) {
    flip <- called[stats::runif(length(called)) < error_rate]
    step <- ifelse(stats::runif(length(flip)) < 0.5, -1L, 1L)
    new_val <- out[flip] + step
    new_val[new_val < 0] <- out[flip][new_val < 0] + 1L
    new_val[new_val > 2] <- out[flip][new_val > 2] - 1L
    out[flip] <- new_val
  }

  if (missing_rate > 0) {
    miss <- called[stats::runif(length(called)) < missing_rate]
    out[miss] <- NA
  }

  out
}
