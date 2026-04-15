<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ACAD-UofA/PolyLinkR.svg?branch=master)](https://travis-ci.org/ACAD-UofA/PolyLinkR)
[![CodeFactor](https://www.codefactor.io/repository/github/acad-uofa/polylinkr/badge)](https://www.codefactor.io/repository/github/acad-uofa/polylinkr)
<!-- badges: end -->

# PolyLink: gene-based pathway enrichment <img src="inst/sticker/polylinkr_150px.png" align="right" />

## Overview

PolyLinkR is an R package that performs gene-based pathway enrichment, which can also be used as evidence for polygenic selection in case the software is used with selection signals evidence. The package explicitly also accounts for linkage disequilibrium between adjacent loci belonging on the same pathway.

PolyLinkR introduces several improvements, and faster implementations of the popular polygenic selection tool [PolySel](https://github.com/CMPG/polysel), which builds upon the core file types and summary statistics used in PolySel. The key difference between PolyLinkR and PolySel is the algorithm used to generate the null distribution of pathway scores. PolySel performs a standard permutation to remap gene scores to genes, whereas PolyLinkR uses a permutation algorithm that randomly links all chromosomes/contigs into a single "circular" genome, and then rotates this circular genome to create a unique mapping between the genes and gene scores. Importantly, this randomisation process preserves the innate linkage structure amongst the genes across the genome, limiting the number of potential false positives that might arise otherwise.

The package exposes a modular pipeline: `plR_read()`, `plR_permute()`, `plR_rescale()`, and `plR_prune()`. Start with `?plR_read` for required inputs (object, set, and mapping tables).

## Installation

Install the development version from GitHub:

```
# install.packages("devtools")
devtools::install_github("ACAD-UofA/PolyLinkR")
```

## Usage

To run pathway or gene set enrichment you need gene-level scores with genomic coordinates, plus tables defining gene sets and gene-to-set membership. The `plR_read()` step assembles these into an internal `plR` object for downstream permutation and rescaling; see the function help pages for argument details and the bundled `data/` objects for small examples.

```r
library(polylinkR)
# ?plR_read
```
