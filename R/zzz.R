# Column names used inside ggplot2::aes() in plot.xpar_boot() are resolved
# via non-standard evaluation against the data frame at plot time, not as
# free variables at package-build time — but static analysis (R CMD check,
# lintr) can't tell that apart from an actual undefined global. Declaring
# them here silences that false positive without adding a hard dependency
# on ggplot2/rlang's `.data` pronoun.
utils::globalVariables(c("rate", "pair_label", "ci_low", "ci_high", "offspring"))
