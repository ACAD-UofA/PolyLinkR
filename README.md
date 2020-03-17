# PolyLinkR
## Overview
PolyLinkR is an R package that performs gene-based pathway enrichment, which can also be used as evidence for polygenic selection in case the software is used with selection signals evidence. The package explicitly also accounts for linkage desiquilibrium between adjacent loci belonging on the same pathway.

PolyLinkR is introduces several improvements, and faster implementations of the popular polygenic selection tool PolySel (github.com/CMPG/polysel), which builds upon the core file types and summary statistics used in PolySel. The key difference between PolyLinkR and PolySel is the algorithm used to generate the null distribution of pathway scores. PolySel performs a standard permutation to remap gene scores to genes, whereas PolyLinkR uses a permutation algorithm that randomly links all chromosomes/contigs into a single 'circular' genome, and then rotates this circular genome to create a unique mapping between the genes and gene scores. Importantly, this randomisation process preserves the innate linkage structure amongst the genes across the genome, limiting the number of potential false positives that might arise otherwise.

## Installation

```
Install from GitHub:
# install.packages("devtools")
devtools::install_github("r-lib/testthat")
```

## Usage



```
output = polylinkr(obj.info = Anatolia_EF_CLR, set.info = PolyLinkR_SetInfo, set.obj = PolyLinkR_SetObj,
             n.cores = 8, emp.nruns = 1000, NN = 10000)
```