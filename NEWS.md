# polylinkR 0.1.0

## Initial CRAN Release

This is the first release of polylinkR to CRAN.

### New Features

* **Gene-based pathway enrichment** with linkage-aware permutation nulls
* Four-step workflow: `plR_read()`, `plR_permute()`, `plR_rescale()`, `plR_prune()`
* **Linkage-preserving permutation** algorithm that maintains chromosome structure
* Support for gene score deconfounding and autocorrelation analysis
* Parallel processing support via `future` framework

### Package Infrastructure

* Comprehensive test suite with 72 tests
* Runnable vignette with example workflow
* Full documentation for all exported functions
* Proper NAMESPACE imports and globalVariables declarations

### Notes

The core statistical algorithms are frozen and validated against the reference implementation in `new_code/`. See the Getting Started vignette for usage examples.
