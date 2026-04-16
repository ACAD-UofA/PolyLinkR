# polylinkR v0.5.0 to v0.6.0 Migration Guide

## Overview

polylinkR v0.6.0 introduces a major naming convention update to align with tidyverse/R package naming standards. This migration guide helps you transition from the legacy naming style (v0.5.0) to the new snake_case convention (v0.6.0).

**Key Changes:**
- Function names: `plR_*` → `*_polylinkr_data` (e.g., `plR_read` → `read_polylinkr_data`)
- Parameter names: dot notation (e.g., `input.path`) → snake_case (e.g., `input_path`)
- Data object names: camelCase (e.g., `obj.info`) → snake_case (e.g., `obj_info`)

---

## Breaking Changes Summary

### 1. Function Name Changes

All main workflow functions have been renamed to follow tidyverse naming conventions:

| Old Name (v0.5.0) | New Name (v0.6.0) | Purpose |
|-------------------|-------------------|---------|
| `plR_read()` | `read_polylinkr_data()` | Read and validate input files |
| `plR_permute()` | `permute_polylinkr_data()` | Perform permutation-based enrichment |
| `plR_rescale()` | `rescale_polylinkr_data()` | Account for genetic autocorrelation |
| `plR_prune()` | `prune_polylinkr_data()` | Adjust for shared genes and FDR correction |

### 2. Parameter Name Changes

#### `read_polylinkr_data()` Parameters

| Old Parameter | New Parameter | Description |
|---------------|---------------|-------------|
| `input.path` | `input_path` | Path to directory with input files |
| `obj.info.path` | `object_info_path` | Path to object info file |
| `set.info.path` | `gene_set_info_path` | Path to gene set info file |
| `set.obj.path` | `gene_set_mapping_path` | Path to gene set mapping file |
| `var.info.path` | `variant_info_path` | Path to variant info file |
| `rec.rate.path` | `recombination_rate_path` | Path to recombination rate file |
| `min.set.n` | `min_set_size` | Minimum gene set size |
| `max.set.n` | `max_set_size` | Maximum gene set size |
| `group` | `group_label` | Label to identify input files |
| `map.fun` | `mapping_function` | Genetic distance mapping function |
| `obj.buffer` | `object_buffer` | Interval around genes (bp) |
| `obj.stat.fun` | `object_statistic_function` | Gene score correction function |
| `bin.size` | `bin_size` | Gene set size interval for correction |
| `obj.in` | `objects_to_include` | Object IDs to retain |
| `obj.out` | `objects_to_exclude` | Object IDs to remove |
| `set.in` | `sets_to_include` | Set IDs to retain |
| `set.out` | `sets_to_exclude` | Set IDs to remove |
| `set.merge` | `merge_threshold` | Minimum proportion for merging |
| `rem.genes` | `remove_duplicate_genes` | Remove identical position genes |

#### `permute_polylinkr_data()` Parameters

| Old Parameter | New Parameter | Description |
|---------------|---------------|-------------|
| `plR.input` | `plr_input` | plR object from previous step |
| `n.perm` | `n_permutations` | Number of permutations |
| `n.boot` | `n_bootstraps` | Number of bootstrap replicates |
| `alt` | `alternative` | Hypothesis test direction |
| `md.meth` | `md_method` | Mahalanobis distance method |
| `kern.wt.max` | `kernel_weight_max` | Maximum probability weight |
| `kern.bound` | `kernel_boundary` | Flanking region for weights |
| `kern.func` | `kernel_function` | Kernel function type |
| `kern.scale` | `kernel_scale` | Kernel scale parameter |
| `gpd.cutoff` | `gpd_cutoff` | GPD tail fitting threshold |
| `n.cores` | `n_cores` | Number of parallel cores |
| `fut.plan` | `future_plan` | Parallel backend strategy |

#### `rescale_polylinkr_data()` Parameters

| Old Parameter | New Parameter | Description |
|---------------|---------------|-------------|
| `plR.input` | `plr_input` | plR object from previous step |
| `cgm.bin` | `cgm_bin` | Minimum pairs per covariance bin |
| `cgm.range` | `cgm_range` | Maximum inter-gene lag |
| `cgm.wt.max` | `cgm_weight_max` | Maximum lag window weight |
| `emp.bayes` | `empirical_bayes` | Empirical Bayes framework |
| `min.rho` | `min_rho` | Minimum correlation threshold |
| `n.cores` | `n_cores` | Number of parallel cores |
| `fut.plan` | `future_plan` | Parallel backend strategy |

#### `prune_polylinkr_data()` Parameters

| Old Parameter | New Parameter | Description |
|---------------|---------------|-------------|
| `plR.input` | `plr_input` | plR object from previous step |
| `n.fdr` | `n_fdr` | Number of FDR replications |
| `est.pi0` | `estimate_pi0` | Estimate pi0 during FDR correction |
| `n.cores` | `n_cores` | Number of parallel cores |
| `fut.plan` | `future_plan` | Parallel backend strategy |

### 3. Data Object Column Name Changes

| Old Name | New Name | Description |
|----------|----------|-------------|
| `objID` | `obj_id` | Unique object (gene) identifier |
| `setID` | `set_id` | Unique gene set identifier |
| `objStat` | `obj_stat` | Object statistic (gene score) |
| `setScore.std` | `set_score_std` | Standardized gene set score |
| `setScore.rs` | `set_score_rs` | Rescaled gene set score |
| `setScore.pr` | `set_score_pr` | Pruned gene set score |
| `setN` | `set_n` | Number of genes in set |

---

## Code Examples: Before vs After

### Example 1: Basic Workflow

#### v0.5.0 (Old)
```r
library(polylinkR)

# Read data
my_plr <- plR_read(
  input.path = "path/to/input",
  min.set.n = 5,
  max.set.n = 500,
  verbose = TRUE
)

# Perform permutation testing
my_plr <- plR_permute(
  plR.input = my_plr,
  n.perm = 1e5,
  n.boot = 30,
  verbose = TRUE
)

# Prune for shared genes
my_plr <- plR_prune(
  plR.input = my_plr,
  n.fdr = 300,
  verbose = TRUE
)
```

#### v0.6.0 (New)
```r
library(polylinkR)

# Read data
my_plr <- read_polylinkr_data(
  input_path = "path/to/input",
  min_set_size = 5,
  max_set_size = 500,
  verbose = TRUE
)

# Perform permutation testing
my_plr <- permute_polylinkr_data(
  plr_input = my_plr,
  n_permutations = 1e5,
  n_bootstraps = 30,
  verbose = TRUE
)

# Prune for shared genes
my_plr <- prune_polylinkr_data(
  plr_input = my_plr,
  n_fdr = 300,
  verbose = TRUE
)
```

### Example 2: Specifying Individual File Paths

#### v0.5.0 (Old)
```r
my_plr <- plR_read(
  obj.info.path = "data/obj.info.csv",
  set.info.path = "data/set.info.csv",
  set.obj.path = "data/set.obj.csv",
  var.info.path = "data/var.info.csv",
  rec.rate.path = "data/rec.rate.csv",
  min.set.n = 2,
  map.fun = "kosambi",
  obj.buffer = 1e4,
  obj.stat.fun = "non.param",
  merge_threshold = 0.95
)
```

#### v0.6.0 (New)
```r
my_plr <- read_polylinkr_data(
  object_info_path = "data/obj.info.csv",
  gene_set_info_path = "data/set.info.csv",
  gene_set_mapping_path = "data/set.obj.csv",
  variant_info_path = "data/var.info.csv",
  recombination_rate_path = "data/rec.rate.csv",
  min_set_size = 2,
  mapping_function = "kosambi",
  object_buffer = 1e4,
  object_statistic_function = "non.param",
  merge_threshold = 0.95
)
```

### Example 3: Advanced Permutation Settings

#### v0.5.0 (Old)
```r
my_plr <- plR_permute(
  plR.input = my_plr,
  n.perm = 5e5,
  n.boot = 50,
  alt = "upper",
  md.meth = "robust",
  kern.wt.max = 0.05,
  kern.bound = 1e5,
  kern.func = "normal",
  kern.scale = 2,
  gpd.cutoff = 0.01,
  n.cores = 4,
  fut.plan = "multisession"
)
```

#### v0.6.0 (New)
```r
my_plr <- permute_polylinkr_data(
  plr_input = my_plr,
  n_permutations = 5e5,
  n_bootstraps = 50,
  alternative = "upper",
  md_method = "robust",
  kernel_weight_max = 0.05,
  kernel_boundary = 1e5,
  kernel_function = "normal",
  kernel_scale = 2,
  gpd_cutoff = 0.01,
  n_cores = 4,
  future_plan = "multisession"
)
```

### Example 4: Rescaling with Custom Parameters

#### v0.5.0 (Old)
```r
my_plr <- plR_rescale(
  plR.input = my_plr,
  rescale = TRUE,
  fast = TRUE,
  cgm.bin = 30,
  cgm.range = 2e6,
  cgm.wt.max = 0.05,
  emp.bayes = "auto",
  min.rho = 1e-5,
  n.cores = 4,
  fut.plan = "multisession"
)
```

#### v0.6.0 (New)
```r
my_plr <- rescale_polylinkr_data(
  plr_input = my_plr,
  rescale = TRUE,
  fast = TRUE,
  cgm_bin = 30,
  cgm_range = 2e6,
  cgm_weight_max = 0.05,
  empirical_bayes = "auto",
  min_rho = 1e-5,
  n_cores = 4,
  future_plan = "multisession"
)
```

---

## Step-by-Step Migration Instructions

### Step 1: Update Package Installation

Ensure you have v0.6.0 or later installed:

```r
# Check current version
packageVersion("polylinkR")

# Install latest version
install.packages("polylinkR")

# Or install from GitHub for development version
# remotes::install_github("ACAD-UofA/PolyLinkR")
```

### Step 2: Identify Legacy Code

Search your codebase for old function and parameter names:

```r
# List of old function names to search for:
# plR_read, plR_permute, plR_rescale, plR_prune

# Common old parameter patterns:
# \.path$, \.n$, \.fun$, \.max$, \.min$, ^obj\., ^set\., ^kern\., ^n\., ^fut\.
```

### Step 3: Update Function Calls

Replace old function names with new ones:

| Old | New |
|-----|-----|
| `plR_read(...)` | `read_polylinkr_data(...)` |
| `plR_permute(...)` | `permute_polylinkr_data(...)` |
| `plR_rescale(...)` | `rescale_polylinkr_data(...)` |
| `plR_prune(...)` | `prune_polylinkr_data(...)` |

### Step 4: Update Parameter Names

Use the mapping tables above to replace parameter names. Common patterns:

1. **Dots to underscores**: `input.path` → `input_path`
2. **Abbreviations expanded**: `n.perm` → `n_permutations`, `alt` → `alternative`
3. **Prefix standardization**: `kern.func` → `kernel_function`, `fut.plan` → `future_plan`
4. **Object prefixes**: `obj.in` → `objects_to_include`, `set.out` → `sets_to_exclude`

### Step 5: Update Data Object References

If you access columns directly from plR objects:

```r
# Old (v0.5.0)
my_plr$obj.info$objID
my_plr$set.info$setScore.std

# New (v0.6.0) - Note: internal data structure names remain as obj.info, set.info
# but column names within may change in future versions
my_plr$obj.info$obj_id
my_plr$set.info$set_score_std
```

### Step 6: Test Your Updated Code

Run your workflow with the new names and verify results:

```r
# Enable warnings to catch any deprecated function usage
options(warn = 1)

# Run your workflow
# ...

# Check for deprecation warnings
```

---

## Suppressing Deprecation Warnings (Temporary)

**Note:** This is for temporary use during migration only. Old names will be removed in v1.0.0.

### Option 1: Suppress All Warnings

```r
# Not recommended - may hide important warnings
options(warn = -1)
# Your code with old names here
options(warn = 0)
```

### Option 2: Use Suppress Warnings (Per Call)

```r
# Suppress deprecation warnings for specific calls
suppressWarnings({
  my_plr <- plR_read(input.path = "data/")
})
```

### Option 3: Capture and Filter Warnings

```r
# Capture warnings and filter out deprecation messages
withCallingHandlers({
  my_plr <- plR_read(input.path = "data/")
}, warning = function(w) {
  if (grepl("deprecated", conditionMessage(w))) {
    invokeRestart("muffleWarning")
  }
})
```

### Option 4: Environment Variable (Advanced)

```r
# Set this before loading polylinkR to suppress deprecation messages
# (Note: This may not be available in all versions)
options(polylinkR.show_deprecation_warnings = FALSE)
library(polylinkR)
```

---

## Timeline for Deprecation

| Version | Status | Legacy Name Support |
|---------|--------|---------------------|
| **v0.5.0** | Previous stable | Full support |
| **v0.6.0** | Current | Deprecated with warnings |
| **v0.7.0** | Future | Deprecated (warnings) |
| **v0.8.0** | Future | Deprecated (warnings) |
| **v0.9.0** | Future | Soft deprecation (reduced warnings) |
| **v1.0.0** | Future | **REMOVED** - Only new names work |

### Recommended Timeline for Users

- **Now (v0.6.0)**: Start planning migration; test new names alongside old
- **Within 3 months**: Complete migration of active projects
- **Within 6 months**: Update all shared/reusable code
- **Before v1.0.0**: Ensure all code uses new names

---

## Troubleshooting

### Issue: "could not find function 'plR_read'"

**Cause**: You're using v0.6.0+ but calling the old function name.

**Solution**: Replace `plR_read` with `read_polylinkr_data`.

### Issue: "unused argument" warnings

**Cause**: You're using old parameter names with new functions.

**Solution**: Check parameter name mapping tables and update.

### Issue: Deprecation warnings appearing in loops

**Cause**: Old names are being called repeatedly.

**Solution**: Update the names inside the loop or wrap in `suppressWarnings()` temporarily.

### Issue: Package scripts not working

**Cause**: External scripts or vignettes may still reference old names.

**Solution**: Update all references or pin dependency to v0.5.0 temporarily:

```r
# In DESCRIPTION or installation
remotes::install_github("ACAD-UofA/PolyLinkR@v0.5.0")
```

---

## Additional Resources

- **NEWS.md**: See `NEWS.md` in the package for detailed change log
- **Function Help**: Use `?read_polylinkr_data` to see updated documentation
- **Vignettes**: Check package vignettes for updated workflow examples
- **GitHub Issues**: https://github.com/ACAD-UofA/PolyLinkR/issues

---

## Quick Reference Card

### Function Names
```
plR_read()           →  read_polylinkr_data()
plR_permute()        →  permute_polylinkr_data()
plR_rescale()        →  rescale_polylinkr_data()
plR_prune()          →  prune_polylinkr_data()
```

### Common Parameter Transformations
```
*.path      →  *_path
*.n         →  *_size / *_count
*.fun       →  *_function
obj.*       →  object_*
set.*       →  set_* / gene_set_*
kern.*      →  kernel_*
n.*         →  n_*
alt         →  alternative
md.meth     →  md_method
fut.plan    →  future_plan
```

---

## Feedback and Support

If you encounter issues during migration:

1. Check this guide for your specific use case
2. Review the updated function documentation with `?function_name`
3. Search existing GitHub issues
4. Open a new issue with:
   - Your code snippet (old and new)
   - Expected vs actual behavior
   - Package version (`packageVersion("polylinkR")`)

---

*Last updated: 2026-04-16*
*For polylinkR v0.6.0*
