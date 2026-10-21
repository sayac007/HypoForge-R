# Validation status and design differences

Native R feature generation, evaluation, safety screening, exploratory ranking,
redundancy screening and rolling-origin Ridge tournament are implemented.

**Not executed in this environment:** R itself and `R CMD check` (R binary unavailable).
No claim of fully verified R functionality or CRAN readiness is made.

The Python ranker uses mutual information; the native R ranker uses absolute
Pearson correlation. The native R selector screens correlation but not the
Python selector's standardized condition number. The R tournament uses native
Ridge normal equations rather than sklearn Ridge, so numerical details may
differ. These are explicitly documented non-parities.

Suggested release gate:
1. `R CMD check .` (zero errors/warnings).
2. `Rscript -e 'testthat::test_local(".")'`.
3. Cross-language golden test vectors for every transformation.
4. Time-series benchmark with an untouched final period.
5. Test installation from clean Windows, Linux, and macOS environments.

No code was copied from third-party R packages; implementation uses base R.
