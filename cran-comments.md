## CRAN Package Submission Comments for polylinkR 0.5.0

### R CMD check results

```
0 errors | 0 warnings | 1 note
```

The single NOTE is:

```
* checking for future file timestamps ... NOTE
  unable to verify current time
```

This is a known issue with CI environments and cannot be fixed. It does not affect package functionality.

### Test Environments

- Local macOS (R 4.4.0): 0 errors, 0 warnings, 1 note
- GitHub Actions Ubuntu (R release): 0 errors
- GitHub Actions Ubuntu (R devel): 0 errors
- GitHub Actions Ubuntu (R oldrel-1): 0 errors
- GitHub Actions macOS (R release): 0 errors
- GitHub Actions Windows (R release): 0 errors

### Downstream Dependencies

None on CRAN.

### Changes in this Version (0.5.0)

* Major documentation enhancement with 5 comprehensive vignettes
* Removed `tdigest` dependency (not available on CRAN), replaced with base R functions
* Added pkgdown site configuration for https://acad-uofa.github.io/PolyLinkR/
* Added GitHub Actions workflows (R-CMD-check, pkgdown, test-coverage)
* Removed Travis CI configuration
* Updated package metadata and documentation

### Method References

There are no published references describing the methods in this package. The package implements original statistical methodology for gene-based pathway enrichment with linkage-aware permutation nulls, as developed by the Australian Centre for Ancient DNA (ACAD) at the University of Adelaide.

### Special Notes

- The package includes a disclaimer in the Description field that it is under active development and should not be used as the sole basis for production or clinical decisions without independent validation.
- All core statistical algorithms (permutation nulls, GPD fitting, rescaling, pruning) are frozen and match the validated reference implementation.
- Test data files in `tests/testthat/humanpops_extract/` are used for integration testing and are not included in the built package.

### URL Checks

All URLs use HTTPS protocol. The GitHub URLs resolve correctly.
