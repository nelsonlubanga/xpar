# xpar

Cross-species parentage assignment for interspecific hybrids.

`xpar` tests which two individuals in a named candidate pool are the parents
of a given offspring, using joint two-parent Mendelian trio consistency.
Unlike tools built around a single, homogeneous population (e.g. `sequoia`,
`CERVUS`, `apparent`), `xpar` does not assume a fixed genotyping-error rate
or Hardy-Weinberg equilibrium. Instead it **self-calibrates**: the expected
discordance rate is measured empirically from a trusted trio (a
pedigree-confirmed pair, or the best-supported pair from a blind search),
and every other candidate pair is judged against that data-derived
reference. This is what lets the same method work whether discordance is
dominated by genotyping noise (a single population) or by real
allele-frequency divergence between parent species (an interspecific
cross).

Uncertainty is quantified with a **block bootstrap over linkage groups**
(not individual SNPs), which is required once markers are correlated by
linkage disequilibrium.

## Installation

```r
# from the package source directory
devtools::install("path/to/xpar")
```

## Basic usage

```r
library(xpar)

# geno: samples x SNPs dosage matrix (0/1/2/NA), scaffold: one label per SNP
d <- xpar_data(geno, scaffold)

# blind scan over every candidate pair
scan <- trio_scan(d, candidates = candidate_ids, offspring = "offspring_id")

# uncertainty + win probability via block bootstrap
boot <- block_bootstrap(d, candidates = candidate_ids, offspring = "offspring_id",
                         n_boot = 1000, seed = 1)
rank_pairs(boot, top_n = 5)

# self-calibrate against a trusted trio, then test a candidate pair against it
ref <- calibrate_reference(d, "trusted_parent1", "trusted_parent2", "offspring_id")
test_vs_reference(rank_pairs(boot, top_n = 1), ref)

plot(boot, top_n = 10)
```

## Status

`0.1.0`. Core trio-consistency engine, empirical calibration, block
bootstrap, group-constrained re-ranking, and plotting are implemented and
unit-tested. An `AlphaSimR`-backed simulation/validation module is also
implemented, sweeping interspecific divergence, genotyping error rate,
missingness, and candidate-pool size (see `vignette("simulation-validation",
package = "xpar")`).
