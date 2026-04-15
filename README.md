  <!-- badges: start -->
  [![R-CMD-check](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/R-CMD-check.yaml)
  [![Travis build status](https://app.travis-ci.com/ACAD-UofA/PolyLinkR.svg?branch=master)](https://app.travis-ci.com/ACAD-UofA/PolyLinkR)
  <!-- badges: end -->
  
# PolyLink: gene-based pathway enrichment <img src="inst/sticker/polylinkr_150px.png" align="right" />

PolyLinkR is an R package that performs gene-based pathway enrichment, which can also be used as evidence for polygenic selection in case the software is used with selection signals evidence. The package explicitly also accounts for linkage desiquilibrium between adjacent loci belonging on the same pathway.

PolyLinkR is introduces several improvements, and faster implementations of the popular polygenic selection tool [PolySel](https://github.com/CMPG/polysel), which builds upon the core file types and summary statistics used in PolySel. The key difference between PolyLinkR and PolySel is the algorithm used to generate the null distribution of pathway scores. PolySel performs a standard permutation to remap gene scores to genes, whereas PolyLinkR uses a permutation algorithm that randomly links all chromosomes/contigs into a single 'circular' genome, and then rotates this circular genome to create a unique mapping between the genes and gene scores. Importantly, this randomisation process preserves the innate linkage structure amongst the genes across the genome, limiting the number of potential false positives that might arise otherwise.

## Installation

Install from CRAN:

```r
install.packages("polylinkR")
```

Or install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("ACAD-UofA/PolyLinkR")
```

## Usage

PolyLinkR uses a four-step workflow:

1. **Read data** with `plR_read()`
2. **Generate permutation nulls** with `plR_permute()`
3. **Rescale for autocorrelation** with `plR_rescale()`
4. **Prune and identify significant sets** with `plR_prune()`

For example:

```r
library(polylinkR)

# Step 1: Read your data
plr_obj <- plR_read(input.path = "path/to/your/data")

# Step 2: Generate permutation nulls
plr_perm <- plR_permute(plr_obj, n.perm = 1e5)

# Step 3: Rescale for genetic autocorrelation
plr_rescale <- plR_rescale(plr_perm)

# Step 4: Prune and identify significant gene sets
plr_final <- plR_prune(plr_rescale)
```

See the [Getting Started vignette](https://acad-uofa.github.io/PolyLinkR/articles/polylinkR.html) for a complete example using the included `tiny_polylinkR` dataset.
