# xpar 0.1.0

* Initial release.
* Core trio-consistency engine (`trio_consistency()`, `trio_scan()`).
* Empirical self-calibration against a trusted trio (`calibrate_reference()`,
  `test_vs_reference()`).
* Block bootstrap over linkage groups for uncertainty quantification
  (`block_bootstrap()`, `rank_pairs()`, `plot.xpar_boot()`).
* Group-constrained re-ranking using independent species/population labels
  (`filter_by_group()`).
* `AlphaSimR`-backed simulation and benchmarking module
  (`simulate_hybrid_pedigree()`, `inject_noise()`, `run_trial()`,
  `benchmark_grid()`), with a pre-computed simulation study (`sim_study`)
  and accompanying vignette.
