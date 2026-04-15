<!-- badges: start -->
[![R-CMD-check](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/R-CMD-check.yaml/badge.svg?branch=dev)](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/pkgdown.yaml/badge.svg?branch=dev)](https://github.com/ACAD-UofA/PolyLinkR/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

# polylinkR: gene-based pathway enrichment with linkage structure <img src="https://raw.githubusercontent.com/ACAD-UofA/PolyLinkR/dev/inst/sticker/polylinkr_150px.png" align="right" alt="" />

## Overview

**polylinkR** performs gene-based pathway enrichment whilst accounting for linkage disequilibrium amongst loci. The package uses a permutation scheme that preserves linkage structure when building null distributions for gene-set scores, which limits false positives relative to unstructured permutations.

The public workflow is modular: `plR_read()`, `plR_permute()`, `plR_rescale()`, and `plR_prune()`. Start with `?plR_read` for required tables and column conventions.

**Active development:** this version is under active development. Do not rely on it as the sole basis for production or clinical decisions without independent validation.

## Documentation site

Function reference, articles, and changelog: <https://acad-uofa.github.io/PolyLinkR/>

## Installation

Until the package is on CRAN, install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("ACAD-UofA/PolyLinkR")
```

Alternatively:

```r
# install.packages("remotes")
remotes::install_github("ACAD-UofA/PolyLinkR")
```

## Local development

- Install [**rig**](https://github.com/r-lib/rig) (optional) for side-by-side **R release** and **R devel** binaries.
- Use [**renv**](https://rstudio.github.io/renv/): from the package root, run `renv::restore()` after cloning. Pak-accelerated installs are enabled in `.Rprofile`; CI sets `RENV_CONFIG_PAK_ENABLED=true` as well.
- Run `devtools::check()` before opening a pull request.

See [`CONTRIBUTING.md`](https://github.com/ACAD-UofA/PolyLinkR/blob/dev/CONTRIBUTING.md) for branch conventions (**`dev`** is the integration branch), release policy, and CRAN preparation notes.

## Related work

The implementation modernises ideas from legacy **PolyLinkR** / **PolyLink** tooling; see the pkgdown article *Getting started with polylinkR* for context and citations.
