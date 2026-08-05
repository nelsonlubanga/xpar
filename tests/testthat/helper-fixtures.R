# Shared fixture used by trio_scan / block_bootstrap / filter_by_group tests.
#
# P1 x P2 -> Off_clean is constructed to be Mendelian-consistent at every
# one of 12 loci (rate = 0), across 3 scaffolds of 4 loci each. P3 is a
# decoy: paired with P1, it produces errors concentrated in scaffolds
# "chr1" and "chr2" but none in "chr3", so the true pair (P1, P2) always
# has the strict lowest rate (0) in every possible bootstrap resample —
# giving a deterministic win_prob == 1 to check against.
make_scan_fixture <- function() {
  p1  <- c(0, 0, 0, 0,   0, 0, 1, 1,   1, 1, 2, 2)
  p2  <- c(0, 0, 1, 1,   0, 2, 1, 1,   2, 2, 2, 2)
  p3  <- c(2, 2, 2, 2,   2, 2, 2, 2,   2, 2, 2, 2)
  off <- c(0, 0, 0, 0,   0, 1, 1, 1,   1, 1, 2, 2)
  scaffold <- rep(c("chr1", "chr2", "chr3"), each = 4)

  geno <- rbind(P1 = p1, P2 = p2, P3 = p3, Off = off)
  colnames(geno) <- paste0("snp", seq_len(ncol(geno)))
  storage.mode(geno) <- "numeric"

  xpar_data(geno, scaffold)
}
