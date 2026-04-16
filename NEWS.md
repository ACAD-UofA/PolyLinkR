# polylinkR 0.5.0

## Major Changes

* **Comprehensive documentation enhancement** with pkgdown site
* Created 5 detailed vignettes:
  - Basic workflow tutorial with Anatolia_EF_CLR dataset
  - Controlling for confounders and covariates
  - Using recombination rates from Bherer et al. 2017
  - Parallel processing for performance optimization
  - Input formats and parameters reference
* Added professional pkgdown configuration with Bootstrap 5
* Added hex logo to package

## Technical Improvements

* Removed `tdigest` dependency (not on CRAN) - replaced with base R `stats::quantile()`
* Added GitHub Actions workflows:
  - R-CMD-check with multi-OS matrix
  - pkgdown auto-deployment to GitHub Pages
  - test-coverage with Codecov integration
* Removed Travis CI (superseded by GitHub Actions)
* Fixed pkgdown reference index for all exported functions

## Package Infrastructure

* 72 comprehensive tests covering core functionality
* Runnable vignettes demonstrating workflows
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

Version 0.5.0 represents a significant documentation milestone. The package implements gene-based pathway enrichment with linkage-preserving permutation nulls. All core statistical algorithms are frozen and validated against the reference implementation.

The pkgdown site is available at: https://acad-uofa.github.io/PolyLinkR/
