## CRAN Package Submission Comments

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
- GitHub Actions Ubuntu (R release + renv): 0 errors
- GitHub Actions Ubuntu (R devel + DESCRIPTION deps): 0 errors

### Downstream Dependencies

None on CRAN yet (first release).

### Method References

There are no published references describing the methods in this package. The package implements original statistical methodology for gene-based pathway enrichment with linkage-aware permutation nulls, as developed by the Australian Centre for Ancient DNA (ACAD) at the University of Adelaide.

### Special Notes

- The package includes a disclaimer in the Description field that it is under active development and should not be used as the sole basis for production or clinical decisions without independent validation.
- All core statistical algorithms (permutation nulls, GPD fitting, rescaling, pruning) are frozen and match the validated reference implementation.
- Test data files in `tests/testthat/humanpops_extract/` are used for integration testing and are not included in the built package.

### URL Checks

All URLs use HTTPS protocol. The GitHub URLs will resolve once the package is published.
