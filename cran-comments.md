## CRAN Package Submission Comments for polylinkR 0.6.0

### R CMD check results

```
0 errors | 0 warnings | 0 notes
```

### Test Environments

- Local macOS (R 4.5.0): 0 errors, 0 warnings, 0 notes
- GitHub Actions Ubuntu (R release): 0 errors
- GitHub Actions Ubuntu (R devel): 0 errors
- GitHub Actions Ubuntu (R oldrel-1): 0 errors
- GitHub Actions macOS (R release): 0 errors
- GitHub Actions Windows (R release): 0 errors

### Downstream Dependencies

None on CRAN.

### Changes in this Version (0.6.0) - Major Update

**Naming Convention Overhaul**

This is a major update from v0.5.0 to v0.6.0 that aligns the package with modern R naming conventions (tidyverse style).

**Function Renaming:**
- All functions have been renamed from camelCase/dot notation to follow tidyverse snake_case conventions
- Example: `PolyLinkR()` → `polylinkr()`, `RunFullPathway()` → `run_full_pathway()`

**Parameter Renaming:**
- All function parameters have been renamed from camelCase/dot notation to snake_case
- Example: `genotypeFile` → `genotype_file`, `pathway.name` → `pathway_name`

**Backward Compatibility:**
- Full backward compatibility is maintained with deprecation warnings
- Old function names remain available but emit `lifecycle::deprecate_warn()` messages
- Old parameter names are still accepted via argument deprecation helpers
- Users are guided to migrate to new naming conventions

**Quality Assurance:**
- All tests pass with the new naming conventions
- All documentation is updated (help files, vignettes, README)
- Examples in all functions use the new snake_case conventions
- No breaking changes for users who follow the migration guidance

### Changes in Previous Version (0.5.0)

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
