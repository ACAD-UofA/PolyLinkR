# polylinkR 0.6.0

## Breaking Changes

This release introduces significant naming convention changes to improve consistency and follow R community standards. **All old names remain functional but emit deprecation warnings** - existing code will continue to work while encouraging migration.

### Function Renaming

* `plR_read()` → `read_polylinkr_data()` - Read and prepare GWAS data
* `plR_permute()` → `permute_polylinkr_data()` - Generate permutation nulls
* `plR_rescale()` → `rescale_polylinkr_data()` - Rescale test statistics
* `plR_prune()` → `prune_polylinkr_data()` - Prune SNPs by linkage

### Parameter Renaming (snake_case convention)

All parameters now use consistent snake_case naming:

| Old Name | New Name |
|----------|----------|
| `plR.input` | `plr_input` |
| `n.perm` | `n_permutations` |
| `n.cores` | `n_cores` |
| `kern.wt.max` | `kernel_weight_max` |
| `kern.wt.type` | `kernel_weight_type` |
| `use.chisq` | `use_chisq` |
| `keep.loci` | `keep_loci` |
| `ld.cutoff` | `ld_cutoff` |
| `ld.info` | `ld_info` |
| `rec.frac.max` | `recombination_fraction_max` |
| `rec.offset` | `recombination_offset` |
| `rec.cM.cut` | `recombination_cutoff_cm` |
| `rec.info` | `recombination_info` |
| `snp2gene.anno` | `snp_to_gene_annotation` |
| `gene2set.anno` | `gene_to_set_annotation` |
| `gene2weight.anno` | `gene_to_weight_annotation` |
| `kernel.fn` | `kernel_function` |
| `seed.val` | `seed_value` |
| `min.gene.set.size` | `min_gene_set_size` |
| `max.gene.set.size` | `max_gene_set_size` |
| `sig.level` | `significance_level` |

### Data Object Renaming

* `Anatolia_EF_CLR` → `anatolia_ef_clr` - Example GWAS dataset
* `PolyLinkR_SetInfo` → `polylinkr_set_info` - Gene set annotation data
* `PolyLinkR_Weights` → `polylinkr_weights` - Gene weight annotations
* `PolyLinkR_SNP2Gene` → `polylinkr_snp2gene` - SNP to gene mapping
* `PolyLinkR_Gene2Set` → `polylinkr_gene2set` - Gene to pathway mapping

### Backward Compatibility

All old function names, parameter names, and data object names **continue to work** and automatically map to the new names with a deprecation warning. This provides a smooth migration path without breaking existing scripts.

**Example migration:**

```r
# Old code (still works with warning)
result <- plR_permute(
  plR.input = data,
  n.perm = 1000,
  n.cores = 4
)

# New code (recommended)
result <- permute_polylinkr_data(
  plr_input = data,
  n_permutations = 1000,
  n_cores = 4
)
```

### Deprecation Timeline

Old names are deprecated and will emit warnings. They are scheduled for removal in a future major release (v1.0.0). Users are encouraged to update their code to use the new naming conventions.

---

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
