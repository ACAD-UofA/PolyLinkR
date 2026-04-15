# polylinkR 0.0.0.9000 (Development Version)

## Package Infrastructure

* Initial CRAN submission preparation.
* Added comprehensive test suite with 72 tests covering core functionality.
* Added runnable vignette demonstrating basic workflow.
* Fixed R CMD check issues:
  - Added missing NAMESPACE imports for stats, utils, methods
  - Added globalVariables for data.table NSE column references
  - Fixed data file organization and documentation
  - Removed stale documentation artifacts
* Fixed code issues:
  - Removed false `plot()` method claims from documentation
  - Fixed missing closing parenthesis in plR_permute.R
  - Fixed `.create_seed()` to return integer values

## Notes

This is a development version preparing for initial CRAN release. The package implements gene-based pathway enrichment with linkage-preserving permutation nulls. All core statistical algorithms are frozen and validated against the reference implementation in `new_code/`.
