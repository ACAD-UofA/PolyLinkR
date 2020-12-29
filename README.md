  <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/ACAD-UofA/PolyLinkR.svg?branch=master)](https://travis-ci.org/ACAD-UofA/PolyLinkR)
  <!-- badges: end -->
  <!-- badges: start -->
[![R build status](https://github.com/ACAD-UofA/PolyLinkR/workflows/R-CMD-check/badge.svg)](https://github.com/ACAD-UofA/PolyLinkR/actions)
<!-- badges: end -->
  
# PolyLink: gene-based pathway enrichment <img src="inst/sticker/polylinkr_150px.png" align="right" />
## Overview
PolyLinkR is an R package that performs gene-based pathway enrichment, which can also be used as evidence for polygenic selection in case the software is used with selection signals evidence. The package explicitly also accounts for linkage desiquilibrium between adjacent loci belonging on the same pathway.

PolyLinkR is introduces several improvements, and faster implementations of the popular polygenic selection tool [PolySel](github.com/CMPG/polysel), which builds upon the core file types and summary statistics used in PolySel. The key difference between PolyLinkR and PolySel is the algorithm used to generate the null distribution of pathway scores. PolySel performs a standard permutation to remap gene scores to genes, whereas PolyLinkR uses a permutation algorithm that randomly links all chromosomes/contigs into a single 'circular' genome, and then rotates this circular genome to create a unique mapping between the genes and gene scores. Importantly, this randomisation process preserves the innate linkage structure amongst the genes across the genome, limiting the number of potential false positives that might arise otherwise.

## Installation

```
Install from GitHub:
# install.packages("devtools")
devtools::install_github("ACAD-UofA/PolyLinkR")
```

## Usage
To run PolyLinkR pathway/gene set enrichment you will need a dataframe/data.table with the geneID, gene scores (from whichever test you're using; e.g. Fst or CLR scores), the chromosome/contig number or ID, and the start/end positions of the gene. Also, two additional dataframes are required, the first definging the gene sets/pathways, and the second linking each gene to the corresponding gene set(s)/pathway(s).

For example, you can use the build in examples:

```
output = polylinkr(obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj,
             n.cores = 8, emp.nruns = 10000, NN = 1000)
```