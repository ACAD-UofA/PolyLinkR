# Getting Started with polylinkR

``` r
library(polylinkR)
```

## Overview

**polylinkR** is an R package for gene-based pathway enrichment that
accounts for linkage disequilibrium among loci. This vignette
demonstrates the basic workflow using the included `tiny_polylinkR`
example dataset.

## The Four-Step Workflow

The polylinkR workflow follows four sequential steps:

1.  **plR_read()** - Read and validate input files
2.  **plR_permute()** - Generate permutation nulls and compute p-values
3.  **plR_rescale()** - Rescale for genetic autocorrelation
4.  **plR_prune()** - Prune gene sets and compute FDR

## Step 1: Read Input Data

First, we read the example data included with the package:

``` r
# Locate the example data
data_path <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
data_path
#> [1] "/home/runner/work/_temp/Library/polylinkR/extdata/tiny_polylinkR"

# Read the data
plr_obj <- plR_read(
  input.path = data_path,
  verbose = FALSE
)

# Check the result
class(plr_obj)
#> [1] "plR"
print(plr_obj)
#> plR class object -- output from plR_read and valid input for plR_permute 
#> NOTE:
#> Standard covariate columns not detected in obj.info
#> ** gene score deconfounding requires valid 'wt.mat' object (see ?polylinkR::plr_permute)
#> ┌────────────────────────────────────────────────────────────────────────────────┐
#> │                                    set.info                                    │
#> └────────────────────────────────────────────────────────────────────────────────┘
#>    setID    setName  setN
#>    <int>     <char> <int>
#> 1:     1 test_set_1     2
#> 2:     2 test_set_2     2
#> ┌────────────────────────────────────────────────────────────────────────────────┐
#> │                                    obj.info                                    │
#> └────────────────────────────────────────────────────────────────────────────────┘
#>    objID objStat objName    chr startpos endpos startpos.base endpos.base
#>    <int>   <num>  <char> <char>    <int>  <int>         <int>       <int>
#> 1:     1     0.5   gene1   chr1     1000   2000          1000        2000
#> 2:     2     0.6   gene2   chr1     3000   4000          3000        4000
#> 3:     3     0.7   gene3   chr1     5000   6000          5000        6000
#>    midpos
#>     <num>
#> 1:   1500
#> 2:   3500
#> 3:   5500
#> ┌────────────────────────────────────────────────────────────────────────────────┐
#> │                                    set.obj                                     │
#> └────────────────────────────────────────────────────────────────────────────────┘
#>    setID objID
#>    <int> <int>
#> 1:     1     1
#> 2:     1     2
#> 3:     2     2
#> 4:     2     3
```

The [`plR_read()`](../reference/plR_read.md) function validates the
input files and creates a `plR` class object containing:

- **set.info**: Gene set information
- **obj.info**: Gene (object) information with scores
- **set.obj**: Gene-to-set mappings

## Inspecting the plR Object

You can inspect the plR object using standard methods:

``` r
# Summary of the object
summary(plr_obj)
#> No tests performed. Nothing to summarise
#> 

# Structure of the object
str(plr_obj)
#> List of 3
#>  $ set.info:Classes 'data.table' and 'data.frame':   2 obs. of  3 variables:
#>   ..$ setID  : int [1:2] 1 2
#>   ..$ setName: chr [1:2] "test_set_1" "test_set_2"
#>   ..$ setN   : int [1:2] 2 2
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ obj.info:Classes 'data.table' and 'data.frame':   3 obs. of  9 variables:
#>   ..$ objID        : int [1:3] 1 2 3
#>   ..$ objStat      : num [1:3] 0.5 0.6 0.7
#>   ..$ objName      : chr [1:3] "gene1" "gene2" "gene3"
#>   ..$ chr          : chr [1:3] "chr1" "chr1" "chr1"
#>   ..$ startpos     : int [1:3] 1000 3000 5000
#>   ..$ endpos       : int [1:3] 2000 4000 6000
#>   ..$ startpos.base: int [1:3] 1000 3000 5000
#>   ..$ endpos.base  : int [1:3] 2000 4000 6000
#>   ..$ midpos       : num [1:3] 1500 3500 5500
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ set.obj :Classes 'data.table' and 'data.frame':   4 obs. of  2 variables:
#>   ..$ setID: int [1:4] 1 1 2 2
#>   ..$ objID: int [1:4] 1 2 2 3
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  - attr(*, "plR.data")=List of 1
#>   ..$ read.data:List of 12
#>   .. ..$ n.genes    : int 3
#>   .. ..$ n.sets     : int 2
#>   .. ..$ n.chr      : int 1
#>   .. ..$ n.cov      : num 0
#>   .. ..$ n.set.genes: int 3
#>   .. ..$ file.paths : Named chr [1:3] "/home/runner/work/_temp/Library/polylinkR/extdata/tiny_polylinkR/ObjInfo.txt" "/home/runner/work/_temp/Library/polylinkR/extdata/tiny_polylinkR/SetInfo.txt" "/home/runner/work/_temp/Library/polylinkR/extdata/tiny_polylinkR/SetObj.txt"
#>   .. .. ..- attr(*, "names")= chr [1:3] "obj.info" "set.info" "set.obj"
#>   .. ..$ coord      : chr "bp"
#>   .. ..$ get.objStat: logi FALSE
#>   .. ..$ cov.names  : NULL
#>   .. ..$ cov.info   : logi FALSE
#>   .. ..$ no.share   : logi FALSE
#>   .. ..$ pos.info   : logi TRUE
#>  - attr(*, "plR.args")=List of 1
#>   ..$ read.args:List of 20
#>   .. ..$ input.path   : chr "/home/runner/work/_temp/Library/polylinkR/extdata/tiny_polylinkR"
#>   .. ..$ obj.info.path: NULL
#>   .. ..$ set.info.path: NULL
#>   .. ..$ set.obj.path : NULL
#>   .. ..$ var.info.path: NULL
#>   .. ..$ rec.rate.path: NULL
#>   .. ..$ min.set.n    : int 2
#>   .. ..$ max.set.n    : num Inf
#>   .. ..$ group        : NULL
#>   .. ..$ map.fun      : chr "kosambi"
#>   .. ..$ obj.buffer   : num 0
#>   .. ..$ obj.stat.fun : chr "non.param"
#>   .. ..$ bin.size     : int 250
#>   .. ..$ obj.in       : NULL
#>   .. ..$ obj.out      : NULL
#>   .. ..$ set.in       : NULL
#>   .. ..$ set.out      : NULL
#>   .. ..$ set.merge    : num 0.95
#>   .. ..$ rem.genes    : logi FALSE
#>   .. ..$ verbose      : logi FALSE
#>  - attr(*, "plR.summary")=List of 1
#>   ..$ read.summary:List of 4
#>   .. ..$ removed.ID    : NULL
#>   .. ..$ repeat.ID     : NULL
#>   .. ..$ merged.ID     : NULL
#>   .. ..$ obj.stat.param: NULL
#>  - attr(*, "plR.session")=List of 1
#>   ..$ read.session:List of 2
#>   .. ..$ session :List of 14
#>   .. .. ..$ R.version  :List of 14
#>   .. .. .. ..$ platform      : chr "x86_64-pc-linux-gnu"
#>   .. .. .. ..$ arch          : chr "x86_64"
#>   .. .. .. ..$ os            : chr "linux-gnu"
#>   .. .. .. ..$ system        : chr "x86_64, linux-gnu"
#>   .. .. .. ..$ status        : chr ""
#>   .. .. .. ..$ major         : chr "4"
#>   .. .. .. ..$ minor         : chr "5.3"
#>   .. .. .. ..$ year          : chr "2026"
#>   .. .. .. ..$ month         : chr "03"
#>   .. .. .. ..$ day           : chr "11"
#>   .. .. .. ..$ svn rev       : chr "89597"
#>   .. .. .. ..$ language      : chr "R"
#>   .. .. .. ..$ version.string: chr "R version 4.5.3 (2026-03-11)"
#>   .. .. .. ..$ nickname      : chr "Reassured Reassurer"
#>   .. .. ..$ platform   : chr "x86_64-pc-linux-gnu"
#>   .. .. ..$ locale     : chr "LC_CTYPE=C.UTF-8;LC_NUMERIC=C;LC_TIME=C.UTF-8;LC_COLLATE=C.UTF-8;LC_MONETARY=C.UTF-8;LC_MESSAGES=C.UTF-8;LC_PAP"| __truncated__
#>   .. .. ..$ tzone      : chr "UTC"
#>   .. .. ..$ tzcode_type: chr "system (glibc)"
#>   .. .. ..$ running    : chr "Ubuntu 24.04.4 LTS"
#>   .. .. ..$ RNGkind    : chr [1:3] "Mersenne-Twister" "Inversion" "Rejection"
#>   .. .. ..$ basePkgs   : chr [1:7] "stats" "graphics" "grDevices" "utils" ...
#>   .. .. ..$ otherPkgs  :List of 1
#>   .. .. .. ..$ polylinkR:List of 24
#>   .. .. .. .. ..$ Package                : chr "polylinkR"
#>   .. .. .. .. ..$ Title                  : chr "Polygenic Pathway Enrichment with Linkage Structure"
#>   .. .. .. .. ..$ Version                : chr "0.5.0"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"ACAD\", \"University of Adelaide\", role = c(\"cph\", \"fnd\")),\n    person(\"PolyLinkR autho"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Perform gene-based pathway enrichment analysis while\n    accounting for linkage disequilibrium (LD) amongst lo"| __truncated__
#>   .. .. .. .. ..$ Language               : chr "en-AU"
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ Roxygen                : chr "list(markdown = TRUE)"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.3"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 4.1.0)"
#>   .. .. .. .. ..$ Imports                : chr "cli, data.table, distances, doFuture, dqrng, foreach, future,\ngstat, igraph, ismev, mgcv, methods, progressr, "| __truncated__
#>   .. .. .. .. ..$ Suggests               : chr "knitr, rmarkdown, spelling, testthat (>= 3.0.0), withr"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ VignetteBuilder        : chr "knitr"
#>   .. .. .. .. ..$ URL                    : chr "https://github.com/ACAD-UofA/PolyLinkR"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/ACAD-UofA/PolyLinkR/issues"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "local::."
#>   .. .. .. .. ..$ RemoteType             : chr "local"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2026-04-16 07:02:53 UTC; runner"
#>   .. .. .. .. ..$ Author                 : chr "ACAD University of Adelaide [cph, fnd],\n  PolyLinkR authors [cre, aut]"
#>   .. .. .. .. ..$ Maintainer             : chr "PolyLinkR authors <polylinkr@acad.edu.au>"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.3; ; 2026-04-16 07:02:55 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/polylinkR/Meta/package.rds"
#>   .. .. ..$ loadedOnly :List of 62
#>   .. .. .. ..$ dqrng       :List of 29
#>   .. .. .. .. ..$ Package          : chr "dqrng"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Fast Pseudo Random Number Generators"
#>   .. .. .. .. ..$ Version          : chr "0.4.1"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n    person(\"Ralf\", \"Stubner\", email = \"ralf.stubner@gmail.com\", role = c(\"aut\", \"cre\"), comment ="| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Several fast random number generators are provided as C++\n  header only libraries: The PCG family by O'Neill ("| __truncated__
#>   .. .. .. .. ..$ License          : chr "AGPL-3"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.5.0)"
#>   .. .. .. .. ..$ Imports          : chr "Rcpp (>= 0.12.16)"
#>   .. .. .. .. ..$ LinkingTo        : chr "Rcpp, BH (>= 1.64.0-1), sitmo (>= 2.0.0)"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.1"
#>   .. .. .. .. ..$ Suggests         : chr "BH, testthat, knitr, rmarkdown, mvtnorm (>= 1.2-3), bench,\nsitmo"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr"
#>   .. .. .. .. ..$ URL              : chr "https://daqana.github.io/dqrng/, https://github.com/daqana/dqrng"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/daqana/dqrng/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2024-05-28 11:40:30 UTC; ralf"
#>   .. .. .. .. ..$ Author           : chr "Ralf Stubner [aut, cre] (<https://orcid.org/0009-0009-1908-106X>),\n  daqana GmbH [cph],\n  David Blackman [cph"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Ralf Stubner <ralf.stubner@gmail.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2024-05-28 22:40:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-01-11 08:09:25 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "dqrng"
#>   .. .. .. .. ..$ RemoteRef        : chr "dqrng"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.4.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/dqrng/Meta/package.rds"
#>   .. .. .. ..$ sass        :List of 29
#>   .. .. .. .. ..$ Type                   : chr "Package"
#>   .. .. .. .. ..$ Package                : chr "sass"
#>   .. .. .. .. ..$ Version                : chr "0.4.10"
#>   .. .. .. .. ..$ Title                  : chr "Syntactically Awesome Style Sheets ('Sass')"
#>   .. .. .. .. ..$ Description            : chr "An 'SCSS' compiler, powered by the 'LibSass' library. With this,\n    R developers can use variables, inheritan"| __truncated__
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Joe\", \"Cheng\", , \"joe@rstudio.com\", \"aut\"),\n    person(\"Timothy\", \"Mastny\", , \"ti"| __truncated__
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://rstudio.github.io/sass/, https://github.com/rstudio/sass"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/rstudio/sass/issues"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.2"
#>   .. .. .. .. ..$ SystemRequirements     : chr "GNU make"
#>   .. .. .. .. ..$ Imports                : chr "fs (>= 1.2.4), rlang (>= 0.4.10), htmltools (>= 0.5.1), R6,\nrappdirs"
#>   .. .. .. .. ..$ Suggests               : chr "testthat, knitr, rmarkdown, withr, shiny, curl"
#>   .. .. .. .. ..$ VignetteBuilder        : chr "knitr"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "yes"
#>   .. .. .. .. ..$ Packaged               : chr "2025-04-11 18:34:19 UTC; cpsievert"
#>   .. .. .. .. ..$ Author                 : chr "Joe Cheng [aut],\n  Timothy Mastny [aut],\n  Richard Iannone [aut] (<https://orcid.org/0000-0003-3925-190X>),\n"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Carson Sievert <carson@rstudio.com>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2025-04-11 19:50:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-14 17:03:33 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "sass"
#>   .. .. .. .. ..$ RemoteRef              : chr "sass"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "0.4.10"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/sass/Meta/package.rds"
#>   .. .. .. ..$ future      :List of 31
#>   .. .. .. .. ..$ Package          : chr "future"
#>   .. .. .. .. ..$ Version          : chr "1.70.0"
#>   .. .. .. .. ..$ Title            : chr "Unified Parallel and Distributed Processing in R for Everyone"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.2.0)"
#>   .. .. .. .. ..$ Imports          : chr "digest, globals (>= 0.18.0), listenv (>= 0.8.0), parallel,\nparallelly (>= 1.44.0), tools, utils"
#>   .. .. .. .. ..$ Suggests         : chr "methods, RhpcBLASctl, R.rsp, markdown"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "R.rsp"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Henrik\", \"Bengtsson\",\n                    role = c(\"aut\", \"cre\", \"cph\"),\n                "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "The purpose of this package is to provide a lightweight and\n    unified Future API for sequential and parallel"| __truncated__
#>   .. .. .. .. ..$ License          : chr "LGPL (>= 2.1)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ ByteCompile      : chr "TRUE"
#>   .. .. .. .. ..$ URL              : chr "https://future.futureverse.org,\nhttps://github.com/futureverse/future"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/future/issues"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ Collate          : chr "'000.bquote.R' '000.import.R' '000.re-exports.R'\n'009.deprecation.R' '010.tweakable.R' '010.utils-parallelly.R"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-13 21:01:52 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0002-7579-5165>)"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-14 06:10:29 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-03-15 05:08:29 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "future"
#>   .. .. .. .. ..$ RemoteRef        : chr "future"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.70.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/future/Meta/package.rds"
#>   .. .. .. ..$ robustbase  :List of 31
#>   .. .. .. .. ..$ Package          : chr "robustbase"
#>   .. .. .. .. ..$ Version          : chr "0.99-7"
#>   .. .. .. .. ..$ VersionNote      : chr "Released 0.99-6 on 2025-09-03, 0.99-5 on 2024-11-01,\n0.99-4-1 on 2024-09-24, 0.99-4 on 2024-08-19 to CRAN"
#>   .. .. .. .. ..$ Date             : chr "2026-02-03"
#>   .. .. .. .. ..$ Title            : chr "Basic Robust Statistics"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Martin\",\"Maechler\", role=c(\"aut\",\"cre\"), email=\"maechler@stat.math.ethz.ch\", comment = c(OR"| __truncated__
#>   .. .. .. .. ..$ URL              : chr "https://robustbase.R-forge.R-project.org/,\nhttps://R-forge.R-project.org/R/?group_id=59,\nhttps://R-forge.R-pr"| __truncated__
#>   .. .. .. .. ..$ BugReports       : chr "https://R-forge.R-project.org/tracker/?atid=302&group_id=59"
#>   .. .. .. .. ..$ Description      : chr "\"Essential\" Robust Statistics.\n Tools allowing to analyze data with robust methods.  This includes\n regress"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.5.0)"
#>   .. .. .. .. ..$ Imports          : chr "stats, graphics, utils, methods, DEoptimR"
#>   .. .. .. .. ..$ Suggests         : chr "grid, MASS, lattice, boot, cluster, Matrix, robust,\nfit.models, MPV, xtable, ggplot2, GGally, RColorBrewer,\nr"| __truncated__
#>   .. .. .. .. ..$ SuggestsNote     : chr "mostly only because of vignette graphics and simulation"
#>   .. .. .. .. ..$ Enhances         : chr "robustX, rrcov, matrixStats, quantreg, Hmisc"
#>   .. .. .. .. ..$ EnhancesNote     : chr "linked to in man/*.Rd"
#>   .. .. .. .. ..$ LazyData         : chr "yes"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ Packaged         : chr "2026-02-04 17:43:48 UTC; maechler"
#>   .. .. .. .. ..$ Author           : chr "Martin Maechler [aut, cre] (ORCID:\n    <https://orcid.org/0000-0002-8685-9910>),\n  Peter Rousseeuw [ctb] (Qn "| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Martin Maechler <maechler@stat.math.ethz.ch>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-02-05 06:10:15 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-02-06 04:46:51 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "robustbase"
#>   .. .. .. .. ..$ RemoteRef        : chr "robustbase"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.99-7"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/robustbase/Meta/package.rds"
#>   .. .. .. ..$ distances   :List of 29
#>   .. .. .. .. ..$ Package          : chr "distances"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Tools for Distance Metrics"
#>   .. .. .. .. ..$ Version          : chr "0.1.13"
#>   .. .. .. .. ..$ Date             : chr "2025-11-23"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Fredrik\", \"Savje\", email = \"rpackages@fredriksavje.com\", role = c(\"aut\", \"cre\")))"
#>   .. .. .. .. ..$ Description      : chr "Provides tools for constructing, manipulating and using distance metrics."
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.4.0)"
#>   .. .. .. .. ..$ Imports          : chr "stats"
#>   .. .. .. .. ..$ Suggests         : chr "testthat"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 3)"
#>   .. .. .. .. ..$ LicenseNote      : chr "The distances packages includes the ANN library\n(distributed under the LGPLv2.1 license)."
#>   .. .. .. .. ..$ URL              : chr "https://github.com/fsavje/distances"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/fsavje/distances/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ Packaged         : chr "2025-11-23 22:12:30 UTC; fredriksavje"
#>   .. .. .. .. ..$ Author           : chr "Fredrik Savje [aut, cre]"
#>   .. .. .. .. ..$ Maintainer       : chr "Fredrik Savje <rpackages@fredriksavje.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-11-24 07:30:10 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-11-25 14:31:15 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "distances"
#>   .. .. .. .. ..$ RemoteRef        : chr "distances"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.1.13"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/distances/Meta/package.rds"
#>   .. .. .. ..$ lattice     :List of 23
#>   .. .. .. .. ..$ Package         : chr "lattice"
#>   .. .. .. .. ..$ Version         : chr "0.22-9"
#>   .. .. .. .. ..$ Date            : chr "2026-02-03"
#>   .. .. .. .. ..$ Priority        : chr "recommended"
#>   .. .. .. .. ..$ Title           : chr "Trellis Graphics for R"
#>   .. .. .. .. ..$ Authors@R       : chr "c(person(\"Deepayan\", \"Sarkar\", role = c(\"aut\", \"cre\"),\n\t            email = \"deepayan.sarkar@r-proje"| __truncated__
#>   .. .. .. .. ..$ Description     : chr "A powerful and elegant high-level data visualization\n  system inspired by Trellis graphics, with an emphasis o"| __truncated__
#>   .. .. .. .. ..$ Depends         : chr "R (>= 4.0.0)"
#>   .. .. .. .. ..$ Suggests        : chr "KernSmooth, MASS, latticeExtra, colorspace"
#>   .. .. .. .. ..$ Imports         : chr "grid, grDevices, graphics, stats, utils"
#>   .. .. .. .. ..$ Enhances        : chr "chron, zoo"
#>   .. .. .. .. ..$ LazyLoad        : chr "yes"
#>   .. .. .. .. ..$ LazyData        : chr "yes"
#>   .. .. .. .. ..$ License         : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL             : chr "https://lattice.r-forge.r-project.org/"
#>   .. .. .. .. ..$ BugReports      : chr "https://github.com/deepayan/lattice/issues"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Packaged        : chr "2026-02-08 19:21:28 UTC; deepayan"
#>   .. .. .. .. ..$ Author          : chr "Deepayan Sarkar [aut, cre] (ORCID:\n    <https://orcid.org/0000-0003-4107-1553>),\n  Felix Andrews [ctb],\n  Ke"| __truncated__
#>   .. .. .. .. ..$ Maintainer      : chr "Deepayan Sarkar <deepayan.sarkar@r-project.org>"
#>   .. .. .. .. ..$ Repository      : chr "CRAN"
#>   .. .. .. .. ..$ Date/Publication: chr "2026-02-09 06:10:13 UTC"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:33:02 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/lattice/Meta/package.rds"
#>   .. .. .. ..$ listenv     :List of 28
#>   .. .. .. .. ..$ Package          : chr "listenv"
#>   .. .. .. .. ..$ Version          : chr "0.10.1"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.1.2)"
#>   .. .. .. .. ..$ Suggests         : chr "R.utils, R.rsp, markdown"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "R.rsp"
#>   .. .. .. .. ..$ Title            : chr "Environments Behaving (Almost) as Lists"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Henrik\", \"Bengtsson\", role=c(\"aut\", \"cre\", \"cph\"),\n                                       "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "List environments are environments that have list-like properties.  For instance, the elements of a list enviro"| __truncated__
#>   .. .. .. .. ..$ License          : chr "LGPL (>= 2.1)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ URL              : chr "https://listenv.futureverse.org,\nhttps://github.com/futureverse/listenv"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/listenv/issues"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-10 16:49:27 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph]"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-10 18:00:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-03-11 04:40:59 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "listenv"
#>   .. .. .. .. ..$ RemoteRef        : chr "listenv"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.10.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/listenv/Meta/package.rds"
#>   .. .. .. ..$ digest      :List of 27
#>   .. .. .. .. ..$ Package          : chr "digest"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Dirk\", \"Eddelbuettel\", role = c(\"aut\", \"cre\"), email = \"edd@debian.org\",\n                 "| __truncated__
#>   .. .. .. .. ..$ Version          : chr "0.6.39"
#>   .. .. .. .. ..$ Date             : chr "2025-11-19"
#>   .. .. .. .. ..$ Title            : chr "Create Compact Hash Digests of R Objects"
#>   .. .. .. .. ..$ Description      : chr "Implementation of a function 'digest()' for the creation of hash\n digests of arbitrary R objects (using the 'm"| __truncated__
#>   .. .. .. .. ..$ URL              : chr "https://github.com/eddelbuettel/digest,\nhttps://eddelbuettel.github.io/digest/,\nhttps://dirk.eddelbuettel.com"| __truncated__
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/eddelbuettel/digest/issues"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.3.0)"
#>   .. .. .. .. ..$ Imports          : chr "utils"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ Suggests         : chr "tinytest, simplermarkdown, rbenchmark"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "simplermarkdown"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2025-11-19 11:55:09 UTC; edd"
#>   .. .. .. .. ..$ Author           : chr "Dirk Eddelbuettel [aut, cre] (ORCID:\n    <https://orcid.org/0000-0001-6419-907X>),\n  Antoine Lucas [ctb] (ORC"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Dirk Eddelbuettel <edd@debian.org>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-11-19 13:20:08 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-11-20 14:16:50 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "digest"
#>   .. .. .. .. ..$ RemoteRef        : chr "digest"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.6.39"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/digest/Meta/package.rds"
#>   .. .. .. ..$ magrittr    :List of 29
#>   .. .. .. .. ..$ Type                : chr "Package"
#>   .. .. .. .. ..$ Package             : chr "magrittr"
#>   .. .. .. .. ..$ Title               : chr "A Forward-Pipe Operator for R"
#>   .. .. .. .. ..$ Version             : chr "2.0.5"
#>   .. .. .. .. ..$ Authors@R           : chr "c(\n    person(\"Stefan Milton\", \"Bache\", , \"stefan@stefanbache.dk\", role = c(\"aut\", \"cph\"),\n        "| __truncated__
#>   .. .. .. .. ..$ Description         : chr "Provides a mechanism for chaining commands with a new\n    forward-pipe operator, %>%. This operator will forwa"| __truncated__
#>   .. .. .. .. ..$ License             : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                 : chr "https://magrittr.tidyverse.org,\nhttps://github.com/tidyverse/magrittr"
#>   .. .. .. .. ..$ BugReports          : chr "https://github.com/tidyverse/magrittr/issues"
#>   .. .. .. .. ..$ Depends             : chr "R (>= 3.4.0)"
#>   .. .. .. .. ..$ Suggests            : chr "covr, knitr, rlang, rmarkdown, testthat"
#>   .. .. .. .. ..$ VignetteBuilder     : chr "knitr"
#>   .. .. .. .. ..$ ByteCompile         : chr "Yes"
#>   .. .. .. .. ..$ Config/Needs/website: chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Encoding            : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote         : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation    : chr "yes"
#>   .. .. .. .. ..$ Packaged            : chr "2026-04-03 08:09:48 UTC; lionel"
#>   .. .. .. .. ..$ Author              : chr "Stefan Milton Bache [aut, cph] (Original author and creator of\n    magrittr),\n  Hadley Wickham [aut],\n  Lion"| __truncated__
#>   .. .. .. .. ..$ Maintainer          : chr "Lionel Henry <lionel@posit.co>"
#>   .. .. .. .. ..$ Repository          : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication    : chr "2026-04-04 08:00:02 UTC"
#>   .. .. .. .. ..$ Built               : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-05 04:44:04 UTC; unix"
#>   .. .. .. .. ..$ RemoteType          : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef        : chr "magrittr"
#>   .. .. .. .. ..$ RemoteRef           : chr "magrittr"
#>   .. .. .. .. ..$ RemoteRepos         : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform   : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha           : chr "2.0.5"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/magrittr/Meta/package.rds"
#>   .. .. .. ..$ evaluate    :List of 28
#>   .. .. .. .. ..$ Type                   : chr "Package"
#>   .. .. .. .. ..$ Package                : chr "evaluate"
#>   .. .. .. .. ..$ Title                  : chr "Parsing and Evaluation Tools that Provide More Details than the\nDefault"
#>   .. .. .. .. ..$ Version                : chr "1.0.5"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Hadley\", \"Wickham\", , \"hadley@posit.co\", role = c(\"aut\", \"cre\")),\n    person(\"Yihui"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Parsing and evaluation tools that make it easy to recreate\n    the command line behaviour of R."
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://evaluate.r-lib.org/, https://github.com/r-lib/evaluate"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/r-lib/evaluate/issues"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 3.6.0)"
#>   .. .. .. .. ..$ Suggests               : chr "callr, covr, ggplot2 (>= 3.3.6), lattice, methods, pkgload,\nragg (>= 1.4.0), rlang (>= 1.1.5), knitr, testthat"| __truncated__
#>   .. .. .. .. ..$ Config/Needs/website   : chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.2"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2025-08-27 16:20:56 UTC; hadleywickham"
#>   .. .. .. .. ..$ Author                 : chr "Hadley Wickham [aut, cre],\n  Yihui Xie [aut] (ORCID: <https://orcid.org/0000-0003-0645-5666>),\n  Michael Lawr"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Hadley Wickham <hadley@posit.co>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2025-08-27 16:40:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2025-08-28 08:11:23 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "evaluate"
#>   .. .. .. .. ..$ RemoteRef              : chr "evaluate"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "1.0.5"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/evaluate/Meta/package.rds"
#>   .. .. .. ..$ grid        :List of 12
#>   .. .. .. .. ..$ Package         : chr "grid"
#>   .. .. .. .. ..$ Version         : chr "4.5.3"
#>   .. .. .. .. ..$ Priority        : chr "base"
#>   .. .. .. .. ..$ Title           : chr "The Grid Graphics Package"
#>   .. .. .. .. ..$ Author          : chr "Paul Murrell <paul@stat.auckland.ac.nz>"
#>   .. .. .. .. ..$ Maintainer      : chr "R Core Team <do-use-Contact-address@r-project.org>"
#>   .. .. .. .. ..$ Contact         : chr "R-help mailing list <r-help@r-project.org>"
#>   .. .. .. .. ..$ Description     : chr "A rewrite of the graphics layout capabilities, plus some\n  support for interaction."
#>   .. .. .. .. ..$ Imports         : chr "grDevices, utils"
#>   .. .. .. .. ..$ License         : chr "Part of R 4.5.3"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:32:05 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/grid/Meta/package.rds"
#>   .. .. .. ..$ iterators   :List of 24
#>   .. .. .. .. ..$ Package          : chr "iterators"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Provides Iterator Construct"
#>   .. .. .. .. ..$ Version          : chr "1.0.14"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Folashade\", \"Daniel\", role=\"cre\", email=\"fdaniel@microsoft.com\"),\n             person(\"Revo"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Support for iterators, which allow a programmer to traverse\n    through all the elements of a vector, list, or"| __truncated__
#>   .. .. .. .. ..$ URL              : chr "https://github.com/RevolutionAnalytics/iterators"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 2.5.0), utils"
#>   .. .. .. .. ..$ Suggests         : chr "RUnit, foreach"
#>   .. .. .. .. ..$ License          : chr "Apache License (== 2.0)"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2022-01-16 18:19:31 UTC; folashade"
#>   .. .. .. .. ..$ Author           : chr "Folashade Daniel [cre],\n  Revolution Analytics [aut, cph],\n  Steve Weston [aut]"
#>   .. .. .. .. ..$ Maintainer       : chr "Folashade Daniel <fdaniel@microsoft.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2022-02-05 00:50:08 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2025-04-05 02:59:54 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "iterators"
#>   .. .. .. .. ..$ RemoteRef        : chr "iterators"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.0.14"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/iterators/Meta/package.rds"
#>   .. .. .. ..$ fastmap     :List of 24
#>   .. .. .. .. ..$ Package          : chr "fastmap"
#>   .. .. .. .. ..$ Title            : chr "Fast Data Structures"
#>   .. .. .. .. ..$ Version          : chr "1.2.0"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n    person(\"Winston\", \"Chang\", email = \"winston@posit.co\", role = c(\"aut\", \"cre\")),\n    person(g"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Fast implementation of data structures, including a key-value\n    store, stack, and queue. Environments are co"| __truncated__
#>   .. .. .. .. ..$ License          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.2.3"
#>   .. .. .. .. ..$ Suggests         : chr "testthat (>= 2.1.1)"
#>   .. .. .. .. ..$ URL              : chr "https://r-lib.github.io/fastmap/, https://github.com/r-lib/fastmap"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/r-lib/fastmap/issues"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2024-05-14 17:54:13 UTC; winston"
#>   .. .. .. .. ..$ Author           : chr "Winston Chang [aut, cre],\n  Posit Software, PBC [cph, fnd],\n  Tessil [cph] (hopscotch_map library)"
#>   .. .. .. .. ..$ Maintainer       : chr "Winston Chang <winston@posit.co>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2024-05-15 09:00:07 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 02:57:16 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "fastmap"
#>   .. .. .. .. ..$ RemoteRef        : chr "fastmap"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.2.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/fastmap/Meta/package.rds"
#>   .. .. .. ..$ Matrix      :List of 27
#>   .. .. .. .. ..$ Package         : chr "Matrix"
#>   .. .. .. .. ..$ Version         : chr "1.7-4"
#>   .. .. .. .. ..$ VersionNote     : chr "do also bump src/version.h, inst/include/Matrix/version.h"
#>   .. .. .. .. ..$ Date            : chr "2025-08-27"
#>   .. .. .. .. ..$ Priority        : chr "recommended"
#>   .. .. .. .. ..$ Title           : chr "Sparse and Dense Matrix Classes and Methods"
#>   .. .. .. .. ..$ Description     : chr "A rich hierarchy of sparse and dense matrix classes,\n\tincluding general, symmetric, triangular, and diagonal "| __truncated__
#>   .. .. .. .. ..$ License         : chr "GPL (>= 2) | file LICENCE"
#>   .. .. .. .. ..$ URL             : chr "https://Matrix.R-forge.R-project.org"
#>   .. .. .. .. ..$ BugReports      : chr "https://R-forge.R-project.org/tracker/?atid=294&group_id=61"
#>   .. .. .. .. ..$ Contact         : chr "Matrix-authors@R-project.org"
#>   .. .. .. .. ..$ Authors@R       : chr "\n\tc(person(\"Douglas\", \"Bates\", role = \"aut\",\n\t         comment = c(ORCID = \"0000-0001-8316-9503\")),"| __truncated__
#>   .. .. .. .. ..$ Depends         : chr "R (>= 4.4), methods"
#>   .. .. .. .. ..$ Imports         : chr "grDevices, graphics, grid, lattice, stats, utils"
#>   .. .. .. .. ..$ Suggests        : chr "MASS, datasets, sfsmisc, tools"
#>   .. .. .. .. ..$ Enhances        : chr "SparseM, graph"
#>   .. .. .. .. ..$ LazyData        : chr "no"
#>   .. .. .. .. ..$ LazyDataNote    : chr "not possible, since we use data/*.R and our S4 classes"
#>   .. .. .. .. ..$ BuildResaveData : chr "no"
#>   .. .. .. .. ..$ Encoding        : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Packaged        : chr "2025-08-27 10:17:11 UTC; maechler"
#>   .. .. .. .. ..$ Author          : chr "Douglas Bates [aut] (ORCID: <https://orcid.org/0000-0001-8316-9503>),\n  Martin Maechler [aut, cre] (ORCID:\n  "| __truncated__
#>   .. .. .. .. ..$ Maintainer      : chr "Martin Maechler <mmaechler+Matrix@gmail.com>"
#>   .. .. .. .. ..$ Repository      : chr "CRAN"
#>   .. .. .. .. ..$ Date/Publication: chr "2025-08-28 11:50:02 UTC"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:33:10 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/Matrix/Meta/package.rds"
#>   .. .. .. ..$ foreach     :List of 29
#>   .. .. .. .. ..$ Package          : chr "foreach"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Provides Foreach Looping Construct"
#>   .. .. .. .. ..$ Version          : chr "1.5.2"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Folashade\", \"Daniel\", role=\"cre\", email=\"fdaniel@microsoft.com\"),\n             person(\"Hong"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Support for the foreach looping construct.  Foreach is an\n        idiom that allows for iterating over element"| __truncated__
#>   .. .. .. .. ..$ License          : chr "Apache License (== 2.0)"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/RevolutionAnalytics/foreach"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/RevolutionAnalytics/foreach/issues"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 2.5.0)"
#>   .. .. .. .. ..$ Imports          : chr "codetools, utils, iterators"
#>   .. .. .. .. ..$ Suggests         : chr "randomForest, doMC, doParallel, testthat, knitr, rmarkdown"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.1.1"
#>   .. .. .. .. ..$ Collate          : chr "'callCombine.R' 'foreach.R' 'do.R' 'foreach-ext.R'\n'foreach-pkg.R' 'getDoPar.R' 'getDoSeq.R' 'getsyms.R' 'iter"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2022-01-11 04:21:12 UTC; fdaniel"
#>   .. .. .. .. ..$ Author           : chr "Folashade Daniel [cre],\n  Hong Ooi [ctb],\n  Rich Calaway [ctb],\n  Microsoft [aut, cph],\n  Steve Weston [aut]"
#>   .. .. .. .. ..$ Maintainer       : chr "Folashade Daniel <fdaniel@microsoft.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2022-02-02 09:20:02 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2025-04-05 03:10:48 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "foreach"
#>   .. .. .. .. ..$ RemoteRef        : chr "foreach"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.5.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/foreach/Meta/package.rds"
#>   .. .. .. ..$ jsonlite    :List of 26
#>   .. .. .. .. ..$ Package          : chr "jsonlite"
#>   .. .. .. .. ..$ Version          : chr "2.0.0"
#>   .. .. .. .. ..$ Title            : chr "A Simple and Robust JSON Parser and Generator for R"
#>   .. .. .. .. ..$ License          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ Depends          : chr "methods"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n    person(\"Jeroen\", \"Ooms\", role = c(\"aut\", \"cre\"), email = \"jeroenooms@gmail.com\",\n      comme"| __truncated__
#>   .. .. .. .. ..$ URL              : chr "https://jeroen.r-universe.dev/jsonlite\nhttps://arxiv.org/abs/1403.2805"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/jeroen/jsonlite/issues"
#>   .. .. .. .. ..$ Maintainer       : chr "Jeroen Ooms <jeroenooms@gmail.com>"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr, R.rsp"
#>   .. .. .. .. ..$ Description      : chr "A reasonably fast JSON parser and generator, optimized for statistical \n    data and the web. Offers simple, f"| __truncated__
#>   .. .. .. .. ..$ Suggests         : chr "httr, vctrs, testthat, knitr, rmarkdown, R.rsp, sf"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.2"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2025-03-26 11:36:10 UTC; jeroen"
#>   .. .. .. .. ..$ Author           : chr "Jeroen Ooms [aut, cre] (<https://orcid.org/0000-0002-4035-0289>),\n  Duncan Temple Lang [ctb],\n  Lloyd Hilaiel"| __truncated__
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-03-27 06:40:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 03:11:11 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "jsonlite"
#>   .. .. .. .. ..$ RemoteRef        : chr "jsonlite"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "2.0.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/jsonlite/Meta/package.rds"
#>   .. .. .. ..$ mgcv        :List of 19
#>   .. .. .. .. ..$ Package         : chr "mgcv"
#>   .. .. .. .. ..$ Version         : chr "1.9-4"
#>   .. .. .. .. ..$ Authors@R       : chr "person(given = \"Simon\",\n                    family = \"Wood\",\n                    role = c(\"aut\", \"cre\"| __truncated__
#>   .. .. .. .. ..$ Title           : chr "Mixed GAM Computation Vehicle with Automatic Smoothness\nEstimation"
#>   .. .. .. .. ..$ Description     : chr "Generalized additive (mixed) models, some of their extensions and \n             other generalized ridge regres"| __truncated__
#>   .. .. .. .. ..$ Priority        : chr "recommended"
#>   .. .. .. .. ..$ Depends         : chr "R (>= 4.4.0), nlme (>= 3.1-64)"
#>   .. .. .. .. ..$ Imports         : chr "methods, stats, graphics, Matrix, splines, utils"
#>   .. .. .. .. ..$ Suggests        : chr "parallel, survival, MASS"
#>   .. .. .. .. ..$ LazyLoad        : chr "yes"
#>   .. .. .. .. ..$ ByteCompile     : chr "yes"
#>   .. .. .. .. ..$ License         : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Packaged        : chr "2025-11-06 12:01:15 UTC; sw283"
#>   .. .. .. .. ..$ Author          : chr "Simon Wood [aut, cre]"
#>   .. .. .. .. ..$ Maintainer      : chr "Simon Wood <simon.wood@r-project.org>"
#>   .. .. .. .. ..$ Repository      : chr "CRAN"
#>   .. .. .. .. ..$ Date/Publication: chr "2025-11-07 16:30:02 UTC"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:36:40 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/mgcv/Meta/package.rds"
#>   .. .. .. ..$ codetools   :List of 15
#>   .. .. .. .. ..$ Package         : chr "codetools"
#>   .. .. .. .. ..$ Version         : chr "0.2-20"
#>   .. .. .. .. ..$ Priority        : chr "recommended"
#>   .. .. .. .. ..$ Author          : chr "Luke Tierney <luke-tierney@uiowa.edu>"
#>   .. .. .. .. ..$ Description     : chr "Code analysis tools for R."
#>   .. .. .. .. ..$ Title           : chr "Code Analysis Tools for R"
#>   .. .. .. .. ..$ Depends         : chr "R (>= 2.1)"
#>   .. .. .. .. ..$ Maintainer      : chr "Luke Tierney <luke-tierney@uiowa.edu>"
#>   .. .. .. .. ..$ URL             : chr "https://gitlab.com/luke-tierney/codetools"
#>   .. .. .. .. ..$ License         : chr "GPL"
#>   .. .. .. .. ..$ NeedsCompilation: chr "no"
#>   .. .. .. .. ..$ Packaged        : chr "2024-03-31 18:18:09 UTC; luke"
#>   .. .. .. .. ..$ Repository      : chr "CRAN"
#>   .. .. .. .. ..$ Date/Publication: chr "2024-03-31 20:10:06 UTC"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; ; 2026-03-11 09:36:22 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/codetools/Meta/package.rds"
#>   .. .. .. ..$ textshaping :List of 32
#>   .. .. .. .. ..$ Package                          : chr "textshaping"
#>   .. .. .. .. ..$ Title                            : chr "Bindings to the 'HarfBuzz' and 'Fribidi' Libraries for Text\nShaping"
#>   .. .. .. .. ..$ Version                          : chr "1.0.5"
#>   .. .. .. .. ..$ Authors@R                        : chr "c(\n    person(\"Thomas Lin\", \"Pedersen\", , \"thomas.pedersen@posit.co\", role = c(\"cre\", \"aut\"),\n     "| __truncated__
#>   .. .. .. .. ..$ Description                      : chr "Provides access to the text shaping functionality in the\n    'HarfBuzz' library and the bidirectional algorith"| __truncated__
#>   .. .. .. .. ..$ License                          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                              : chr "https://github.com/r-lib/textshaping"
#>   .. .. .. .. ..$ BugReports                       : chr "https://github.com/r-lib/textshaping/issues"
#>   .. .. .. .. ..$ Depends                          : chr "R (>= 3.2.0)"
#>   .. .. .. .. ..$ Imports                          : chr "lifecycle, stats, stringi, systemfonts (>= 1.3.0), utils"
#>   .. .. .. .. ..$ Suggests                         : chr "covr, grDevices, grid, knitr, rmarkdown, testthat (>= 3.0.0)"
#>   .. .. .. .. ..$ LinkingTo                        : chr "cpp11 (>= 0.2.1), systemfonts (>= 1.0.0)"
#>   .. .. .. .. ..$ VignetteBuilder                  : chr "knitr"
#>   .. .. .. .. ..$ Config/build/compilation-database: chr "true"
#>   .. .. .. .. ..$ Config/testthat/edition          : chr "3"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep       : chr "2025-04-23"
#>   .. .. .. .. ..$ Encoding                         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                      : chr "7.3.2"
#>   .. .. .. .. ..$ SystemRequirements               : chr "freetype2, harfbuzz, fribidi"
#>   .. .. .. .. ..$ NeedsCompilation                 : chr "yes"
#>   .. .. .. .. ..$ Packaged                         : chr "2026-03-06 09:57:58 UTC; thomas"
#>   .. .. .. .. ..$ Author                           : chr "Thomas Lin Pedersen [cre, aut] (ORCID:\n    <https://orcid.org/0000-0002-5147-4711>),\n  Posit Software, PBC [c"| __truncated__
#>   .. .. .. .. ..$ Maintainer                       : chr "Thomas Lin Pedersen <thomas.pedersen@posit.co>"
#>   .. .. .. .. ..$ Repository                       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication                 : chr "2026-03-06 10:40:02 UTC"
#>   .. .. .. .. ..$ Built                            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-05 04:50:20 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                     : chr "textshaping"
#>   .. .. .. .. ..$ RemoteRef                        : chr "textshaping"
#>   .. .. .. .. ..$ RemoteRepos                      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform                : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                        : chr "1.0.5"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/textshaping/Meta/package.rds"
#>   .. .. .. ..$ jquerylib   :List of 24
#>   .. .. .. .. ..$ Package                : chr "jquerylib"
#>   .. .. .. .. ..$ Title                  : chr "Obtain 'jQuery' as an HTML Dependency Object"
#>   .. .. .. .. ..$ Version                : chr "0.1.4"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Carson\", \"Sievert\", role = c(\"aut\", \"cre\"), email = \"carson@rstudio.com\", comment = c"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Obtain any major version of 'jQuery' (<https://code.jquery.com/>) and use it in any webpage generated by 'htmlt"| __truncated__
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.0.2"
#>   .. .. .. .. ..$ Imports                : chr "htmltools"
#>   .. .. .. .. ..$ Suggests               : chr "testthat"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2021-04-26 16:40:21 UTC; cpsievert"
#>   .. .. .. .. ..$ Author                 : chr "Carson Sievert [aut, cre] (<https://orcid.org/0000-0002-4958-2844>),\n  Joe Cheng [aut],\n  RStudio [cph],\n  j"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Carson Sievert <carson@rstudio.com>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2021-04-26 17:10:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2025-04-05 03:28:11 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "jquerylib"
#>   .. .. .. .. ..$ RemoteRef              : chr "jquerylib"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "0.1.4"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/jquerylib/Meta/package.rds"
#>   .. .. .. ..$ cli         :List of 29
#>   .. .. .. .. ..$ Package                   : chr "cli"
#>   .. .. .. .. ..$ Title                     : chr "Helpers for Developing Command Line Interfaces"
#>   .. .. .. .. ..$ Version                   : chr "3.6.6"
#>   .. .. .. .. ..$ Authors@R                 : chr "c(\n    person(\"Gábor\", \"Csárdi\", , \"gabor@posit.co\", role = c(\"aut\", \"cre\")),\n    person(\"Hadley\""| __truncated__
#>   .. .. .. .. ..$ Description               : chr "A suite of tools to build attractive command line interfaces\n    ('CLIs'), from semantic elements: headings, l"| __truncated__
#>   .. .. .. .. ..$ License                   : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                       : chr "https://cli.r-lib.org, https://github.com/r-lib/cli"
#>   .. .. .. .. ..$ BugReports                : chr "https://github.com/r-lib/cli/issues"
#>   .. .. .. .. ..$ Depends                   : chr "R (>= 3.4)"
#>   .. .. .. .. ..$ Imports                   : chr "utils"
#>   .. .. .. .. ..$ Suggests                  : chr "callr, covr, crayon, digest, glue (>= 1.6.0), grDevices,\nhtmltools, htmlwidgets, knitr, methods, processx, ps "| __truncated__
#>   .. .. .. .. ..$ Config/Needs/website      : chr "r-lib/asciicast, bench, brio, cpp11, decor, desc,\nfansi, prettyunits, sessioninfo, tidyverse/tidytemplate,\nusethis, vctrs"
#>   .. .. .. .. ..$ Config/testthat/edition   : chr "3"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep: chr "2025-04-25"
#>   .. .. .. .. ..$ Encoding                  : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote               : chr "7.3.2.9000"
#>   .. .. .. .. ..$ NeedsCompilation          : chr "yes"
#>   .. .. .. .. ..$ Packaged                  : chr "2026-04-08 18:14:45 UTC; gaborcsardi"
#>   .. .. .. .. ..$ Author                    : chr "Gábor Csárdi [aut, cre],\n  Hadley Wickham [ctb],\n  Kirill Müller [ctb],\n  Salim Brüggemann [ctb] (ORCID: <ht"| __truncated__
#>   .. .. .. .. ..$ Maintainer                : chr "Gábor Csárdi <gabor@posit.co>"
#>   .. .. .. .. ..$ Repository                : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication          : chr "2026-04-09 09:50:18 UTC"
#>   .. .. .. .. ..$ Built                     : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-10 04:38:10 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef              : chr "cli"
#>   .. .. .. .. ..$ RemoteRef                 : chr "cli"
#>   .. .. .. .. ..$ RemoteRepos               : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform         : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                 : chr "3.6.6"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/cli/Meta/package.rds"
#>   .. .. .. ..$ rlang       :List of 32
#>   .. .. .. .. ..$ Package                          : chr "rlang"
#>   .. .. .. .. ..$ Version                          : chr "1.2.0"
#>   .. .. .. .. ..$ Title                            : chr "Functions for Base Types and Core R and 'Tidyverse' Features"
#>   .. .. .. .. ..$ Description                      : chr "A toolbox for working with base types, core R features\n  like the condition system, and core 'Tidyverse' featu"| __truncated__
#>   .. .. .. .. ..$ Authors@R                        : chr "c(\n    person(\"Lionel\", \"Henry\", ,\"lionel@posit.co\", c(\"aut\", \"cre\")),\n    person(\"Hadley\", \"Wic"| __truncated__
#>   .. .. .. .. ..$ License                          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ ByteCompile                      : chr "true"
#>   .. .. .. .. ..$ Biarch                           : chr "true"
#>   .. .. .. .. ..$ Depends                          : chr "R (>= 4.0.0)"
#>   .. .. .. .. ..$ Imports                          : chr "utils"
#>   .. .. .. .. ..$ Suggests                         : chr "cli (>= 3.1.0), covr, crayon, desc, fs, glue, knitr,\nmagrittr, methods, pillar, pkgload, rmarkdown, stats, tes"| __truncated__
#>   .. .. .. .. ..$ Enhances                         : chr "winch"
#>   .. .. .. .. ..$ Encoding                         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                      : chr "7.3.3"
#>   .. .. .. .. ..$ URL                              : chr "https://rlang.r-lib.org, https://github.com/r-lib/rlang"
#>   .. .. .. .. ..$ BugReports                       : chr "https://github.com/r-lib/rlang/issues"
#>   .. .. .. .. ..$ Config/build/compilation-database: chr "true"
#>   .. .. .. .. ..$ Config/testthat/edition          : chr "3"
#>   .. .. .. .. ..$ Config/Needs/website             : chr "dplyr, tidyverse/tidytemplate"
#>   .. .. .. .. ..$ NeedsCompilation                 : chr "yes"
#>   .. .. .. .. ..$ Packaged                         : chr "2026-04-02 12:23:10 UTC; lionel"
#>   .. .. .. .. ..$ Author                           : chr "Lionel Henry [aut, cre],\n  Hadley Wickham [aut],\n  mikefc [cph] (Hash implementation based on Mike's xxhashli"| __truncated__
#>   .. .. .. .. ..$ Maintainer                       : chr "Lionel Henry <lionel@posit.co>"
#>   .. .. .. .. ..$ Repository                       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication                 : chr "2026-04-06 10:40:02 UTC"
#>   .. .. .. .. ..$ Built                            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-07 04:36:36 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                     : chr "rlang"
#>   .. .. .. .. ..$ RemoteRef                        : chr "rlang"
#>   .. .. .. .. ..$ RemoteRepos                      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform                : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                        : chr "1.2.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/rlang/Meta/package.rds"
#>   .. .. .. ..$ zigg        :List of 25
#>   .. .. .. .. ..$ Package          : chr "zigg"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Lightweight Interfaces to the 'Ziggurat' Pseudo Random Number\nGenerator"
#>   .. .. .. .. ..$ Version          : chr "0.0.2"
#>   .. .. .. .. ..$ Date             : chr "2025-02-07"
#>   .. .. .. .. ..$ Authors@R        : chr "person(\"Dirk\", \"Eddelbuettel\", role = c(\"aut\", \"cre\"), email = \"edd@debian.org\",\n                  c"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "The 'Ziggurat' pseudo-random number generator (or PRNG),\n introduced by Marsaglia and Tsang (2000, <doi:10.186"| __truncated__
#>   .. .. .. .. ..$ URL              : chr "https://github.com/eddelbuettel/zigg"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/eddelbuettel/zigg/issues"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "6.0.1"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2025-02-07 13:59:06 UTC; edd"
#>   .. .. .. .. ..$ Author           : chr "Dirk Eddelbuettel [aut, cre] (<https://orcid.org/0000-0001-6419-907X>)"
#>   .. .. .. .. ..$ Maintainer       : chr "Dirk Eddelbuettel <edd@debian.org>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-02-07 14:20:01 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 03:17:53 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "zigg"
#>   .. .. .. .. ..$ RemoteRef        : chr "zigg"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.0.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/zigg/Meta/package.rds"
#>   .. .. .. ..$ parallelly  :List of 29
#>   .. .. .. .. ..$ Package          : chr "parallelly"
#>   .. .. .. .. ..$ Version          : chr "1.46.1"
#>   .. .. .. .. ..$ Title            : chr "Enhancing the 'parallel' Package"
#>   .. .. .. .. ..$ Imports          : chr "parallel, tools, utils"
#>   .. .. .. .. ..$ Suggests         : chr "commonmark, base64enc"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "parallelly"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n        person(\"Henrik\", \"Bengtsson\",\n               role = c(\"aut\", \"cre\", \"cph\"),\n           "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Utility functions that enhance the 'parallel' package and support the built-in parallel backends of the 'future"| __truncated__
#>   .. .. .. .. ..$ License          : chr "LGPL (>= 2.1)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ ByteCompile      : chr "TRUE"
#>   .. .. .. .. ..$ URL              : chr "https://parallelly.futureverse.org,\nhttps://github.com/futureverse/parallelly"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/parallelly/issues"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-01-07 07:32:26 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0002-7579-5165>),\n  Mike Cheng [ctb]"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-01-08 06:10:23 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-01-09 04:36:07 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "parallelly"
#>   .. .. .. .. ..$ RemoteRef        : chr "parallelly"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.46.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/parallelly/Meta/package.rds"
#>   .. .. .. ..$ future.apply:List of 29
#>   .. .. .. .. ..$ Package          : chr "future.apply"
#>   .. .. .. .. ..$ Version          : chr "1.20.2"
#>   .. .. .. .. ..$ Title            : chr "Apply Function to Elements in Parallel using Futures"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.2.0), future (>= 1.49.0)"
#>   .. .. .. .. ..$ Imports          : chr "globals, parallel, utils"
#>   .. .. .. .. ..$ Suggests         : chr "datasets, stats, tools, listenv, R.rsp, markdown"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "R.rsp"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Henrik\", \"Bengtsson\",\n                    role = c(\"aut\", \"cre\", \"cph\"),\n                "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Implementations of apply(), by(), eapply(), lapply(), Map(), .mapply(), mapply(), replicate(), sapply(), tapply"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ URL              : chr "https://future.apply.futureverse.org,\nhttps://github.com/futureverse/future.apply"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/future.apply/issues"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-02-20 04:06:25 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0002-7579-5165>),\n  R Core Team [cph, ctb]"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-02-20 12:00:09 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-02-21 04:31:11 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "future.apply"
#>   .. .. .. .. ..$ RemoteRef        : chr "future.apply"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.20.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/future.apply/Meta/package.rds"
#>   .. .. .. ..$ splines     :List of 13
#>   .. .. .. .. ..$ Package         : chr "splines"
#>   .. .. .. .. ..$ Version         : chr "4.5.3"
#>   .. .. .. .. ..$ Priority        : chr "base"
#>   .. .. .. .. ..$ Imports         : chr "graphics, stats"
#>   .. .. .. .. ..$ Title           : chr "Regression Spline Functions and Classes"
#>   .. .. .. .. ..$ Author          : chr "Douglas M. Bates <bates@stat.wisc.edu> and\n William N. Venables <Bill.Venables@csiro.au>"
#>   .. .. .. .. ..$ Maintainer      : chr "R Core Team <do-use-Contact-address@r-project.org>"
#>   .. .. .. .. ..$ Contact         : chr "R-help mailing list <r-help@r-project.org>"
#>   .. .. .. .. ..$ Description     : chr "Regression spline functions and classes."
#>   .. .. .. .. ..$ License         : chr "Part of R 4.5.3"
#>   .. .. .. .. ..$ Suggests        : chr "Matrix, methods"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:32:12 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/splines/Meta/package.rds"
#>   .. .. .. ..$ intervals   :List of 25
#>   .. .. .. .. ..$ Package          : chr "intervals"
#>   .. .. .. .. ..$ Version          : chr "0.15.5"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Tools for Working with Points and Intervals"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(given = \"Richard\",\n                    family = \"Bourgon\",\n                    role = \"aut\",\n"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 2.9.0)"
#>   .. .. .. .. ..$ Imports          : chr "utils, graphics, methods"
#>   .. .. .. .. ..$ Description      : chr "Tools for working with and comparing sets of points and intervals."
#>   .. .. .. .. ..$ License          : chr "Artistic-2.0"
#>   .. .. .. .. ..$ LazyLoad         : chr "yes"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/edzer/intervals"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2024-08-23 10:51:00 UTC; edzer"
#>   .. .. .. .. ..$ Author           : chr "Richard Bourgon [aut],\n  Edzer Pebesma [cre]"
#>   .. .. .. .. ..$ Maintainer       : chr "Edzer Pebesma <edzer.pebesma@uni-muenster.de>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2024-08-23 12:10:01 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 02:59:45 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "intervals"
#>   .. .. .. .. ..$ RemoteRef        : chr "intervals"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.15.5"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/intervals/Meta/package.rds"
#>   .. .. .. ..$ cachem      :List of 27
#>   .. .. .. .. ..$ Package             : chr "cachem"
#>   .. .. .. .. ..$ Version             : chr "1.1.0"
#>   .. .. .. .. ..$ Title               : chr "Cache R Objects with Automatic Pruning"
#>   .. .. .. .. ..$ Description         : chr "Key-value stores with automatic pruning. Caches can limit\n    either their total size or the age of the oldest"| __truncated__
#>   .. .. .. .. ..$ Authors@R           : chr "c(\n    person(\"Winston\", \"Chang\", , \"winston@posit.co\", c(\"aut\", \"cre\")),\n    person(family = \"Pos"| __truncated__
#>   .. .. .. .. ..$ License             : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ Encoding            : chr "UTF-8"
#>   .. .. .. .. ..$ ByteCompile         : chr "true"
#>   .. .. .. .. ..$ URL                 : chr "https://cachem.r-lib.org/, https://github.com/r-lib/cachem"
#>   .. .. .. .. ..$ Imports             : chr "rlang, fastmap (>= 1.2.0)"
#>   .. .. .. .. ..$ Suggests            : chr "testthat"
#>   .. .. .. .. ..$ RoxygenNote         : chr "7.2.3"
#>   .. .. .. .. ..$ Config/Needs/routine: chr "lobstr"
#>   .. .. .. .. ..$ Config/Needs/website: chr "pkgdown"
#>   .. .. .. .. ..$ NeedsCompilation    : chr "yes"
#>   .. .. .. .. ..$ Packaged            : chr "2024-05-15 15:54:22 UTC; winston"
#>   .. .. .. .. ..$ Author              : chr "Winston Chang [aut, cre],\n  Posit Software, PBC [cph, fnd]"
#>   .. .. .. .. ..$ Maintainer          : chr "Winston Chang <winston@posit.co>"
#>   .. .. .. .. ..$ Repository          : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication    : chr "2024-05-16 09:50:11 UTC"
#>   .. .. .. .. ..$ Built               : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 03:12:33 UTC; unix"
#>   .. .. .. .. ..$ RemoteType          : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef        : chr "cachem"
#>   .. .. .. .. ..$ RemoteRef           : chr "cachem"
#>   .. .. .. .. ..$ RemoteRepos         : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform   : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha           : chr "1.1.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/cachem/Meta/package.rds"
#>   .. .. .. ..$ yaml        :List of 28
#>   .. .. .. .. ..$ Type                   : chr "Package"
#>   .. .. .. .. ..$ Package                : chr "yaml"
#>   .. .. .. .. ..$ Title                  : chr "Methods to Convert R Data to YAML and Back"
#>   .. .. .. .. ..$ Version                : chr "2.3.12"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Hadley\", \"Wickham\", , \"hadley@posit.co\", role = \"cre\",\n           comment = c(ORCID = "| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Implements the 'libyaml' 'YAML' 1.1 parser and emitter\n    (<https://pyyaml.org/wiki/LibYAML>) for R."
#>   .. .. .. .. ..$ License                : chr "BSD_3_clause + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://yaml.r-lib.org, https://github.com/r-lib/yaml/"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/r-lib/yaml/issues"
#>   .. .. .. .. ..$ Suggests               : chr "knitr, rmarkdown, testthat (>= 3.0.0)"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Config/Needs/website   : chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.3"
#>   .. .. .. .. ..$ VignetteBuilder        : chr "knitr"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "yes"
#>   .. .. .. .. ..$ Packaged               : chr "2025-12-08 16:53:07 UTC; hadleywickham"
#>   .. .. .. .. ..$ Author                 : chr "Hadley Wickham [cre] (ORCID: <https://orcid.org/0000-0003-4757-117X>),\n  Shawn Garbett [ctb] (ORCID: <https://"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Hadley Wickham <hadley@posit.co>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2025-12-10 07:00:01 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-12-10 19:48:14 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "yaml"
#>   .. .. .. .. ..$ RemoteRef              : chr "yaml"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "2.3.12"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/yaml/Meta/package.rds"
#>   .. .. .. ..$ FNN         :List of 24
#>   .. .. .. .. ..$ Package          : chr "FNN"
#>   .. .. .. .. ..$ Version          : chr "1.1.4.1"
#>   .. .. .. .. ..$ Date             : chr "2023-12-31"
#>   .. .. .. .. ..$ Title            : chr "Fast Nearest Neighbor Search Algorithms and Applications"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Alina\", \"Beygelzimer\", role = \"aut\",\n                    comment = \"cover tree library\"),\n "| __truncated__
#>   .. .. .. .. ..$ Copyright        : chr "ANN Copyright (c) 1997-2010 University of Maryland and Sunil\nArya and David Mount. All Rights Reserved."
#>   .. .. .. .. ..$ Depends          : chr "R (>= 4.0.0)"
#>   .. .. .. .. ..$ Suggests         : chr "chemometrics, mvtnorm"
#>   .. .. .. .. ..$ Description      : chr "Cover-tree and kd-tree fast k-nearest neighbor search algorithms and related applications\n        including KN"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2024-09-22 07:57:30 UTC; hornik"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2024-09-22 08:24:48 UTC"
#>   .. .. .. .. ..$ Author           : chr "Alina Beygelzimer [aut] (cover tree library),\n  Sham Kakadet [aut] (cover tree library),\n  John Langford [aut"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Shengqiao Li <lishengqiao@yahoo.com>"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-04-05 02:56:22 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "FNN"
#>   .. .. .. .. ..$ RemoteRef        : chr "FNN"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.1.4.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/FNN/Meta/package.rds"
#>   .. .. .. ..$ tools       :List of 12
#>   .. .. .. .. ..$ Package         : chr "tools"
#>   .. .. .. .. ..$ Version         : chr "4.5.3"
#>   .. .. .. .. ..$ Priority        : chr "base"
#>   .. .. .. .. ..$ Title           : chr "Tools for Package Development"
#>   .. .. .. .. ..$ Author          : chr "R Core Team"
#>   .. .. .. .. ..$ Maintainer      : chr "R Core Team <do-use-Contact-address@r-project.org>"
#>   .. .. .. .. ..$ Contact         : chr "R-help mailing list <r-help@r-project.org>"
#>   .. .. .. .. ..$ Description     : chr "Tools for package development, administration and documentation."
#>   .. .. .. .. ..$ License         : chr "Part of R 4.5.3"
#>   .. .. .. .. ..$ Suggests        : chr "codetools, methods, xml2, curl, commonmark, knitr, xfun, mathjaxr, V8"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:30:23 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/tools/Meta/package.rds"
#>   .. .. .. ..$ parallel    :List of 14
#>   .. .. .. .. ..$ Package         : chr "parallel"
#>   .. .. .. .. ..$ Version         : chr "4.5.3"
#>   .. .. .. .. ..$ Priority        : chr "base"
#>   .. .. .. .. ..$ Title           : chr "Support for Parallel Computation in R"
#>   .. .. .. .. ..$ Author          : chr "R Core Team"
#>   .. .. .. .. ..$ Maintainer      : chr "R Core Team <do-use-Contact-address@r-project.org>"
#>   .. .. .. .. ..$ Contact         : chr "R-help mailing list <r-help@r-project.org>"
#>   .. .. .. .. ..$ Description     : chr "Support for parallel computation, including by forking\n   (taken from package multicore), by sockets (taken fr"| __truncated__
#>   .. .. .. .. ..$ License         : chr "Part of R 4.5.3"
#>   .. .. .. .. ..$ Imports         : chr "tools, compiler"
#>   .. .. .. .. ..$ Suggests        : chr "methods"
#>   .. .. .. .. ..$ Enhances        : chr "snow, Rmpi, mirai"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:32:15 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/parallel/Meta/package.rds"
#>   .. .. .. ..$ doFuture    :List of 29
#>   .. .. .. .. ..$ Package          : chr "doFuture"
#>   .. .. .. .. ..$ Version          : chr "1.2.1"
#>   .. .. .. .. ..$ Title            : chr "Use Foreach to Parallelize via the Future Framework"
#>   .. .. .. .. ..$ Depends          : chr "foreach (>= 1.5.0), future (>= 1.49.0)"
#>   .. .. .. .. ..$ Imports          : chr "future.apply, globals, iterators, parallel, utils"
#>   .. .. .. .. ..$ Suggests         : chr "doRNG (>= 1.8.2), markdown, R.rsp"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "R.rsp"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Henrik\", \"Bengtsson\",\n                    role = c(\"aut\", \"cre\", \"cph\"),\n                "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "The 'future' package provides a unifying parallelization framework for R that supports many parallel and distri"| __truncated__
#>   .. .. .. .. ..$ License          : chr "LGPL (>= 2.1)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ URL              : chr "https://doFuture.futureverse.org,\nhttps://github.com/futureverse/doFuture"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/doFuture/issues"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-02-19 22:55:07 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0002-7579-5165>)"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-02-20 06:40:17 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-02-21 04:31:36 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "doFuture"
#>   .. .. .. .. ..$ RemoteRef        : chr "doFuture"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.2.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/doFuture/Meta/package.rds"
#>   .. .. .. ..$ globals     :List of 28
#>   .. .. .. .. ..$ Package          : chr "globals"
#>   .. .. .. .. ..$ Version          : chr "0.19.1"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.1.2)"
#>   .. .. .. .. ..$ Imports          : chr "codetools"
#>   .. .. .. .. ..$ Title            : chr "Identify Global Objects in R Expressions"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n    person(\"Henrik\", \"Bengtsson\", role=c(\"aut\", \"cre\", \"cph\"),\n           email=\"henrikb@braju."| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Identifies global (\"unknown\" or \"free\") objects in R expressions\n    by code inspection using various stra"| __truncated__
#>   .. .. .. .. ..$ License          : chr "LGPL (>= 2.1)"
#>   .. .. .. .. ..$ LazyLoad         : chr "TRUE"
#>   .. .. .. .. ..$ ByteCompile      : chr "TRUE"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ URL              : chr "https://globals.futureverse.org,\nhttps://github.com/futureverse/globals"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/globals/issues"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-12 23:06:06 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph],\n  Davis Vaughan [ctb]"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-13 07:30:09 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-03-14 04:31:23 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "globals"
#>   .. .. .. .. ..$ RemoteRef        : chr "globals"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.19.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/globals/Meta/package.rds"
#>   .. .. .. ..$ Rfast       :List of 28
#>   .. .. .. .. ..$ Package           : chr "Rfast"
#>   .. .. .. .. ..$ Type              : chr "Package"
#>   .. .. .. .. ..$ Title             : chr "A Collection of Efficient and Extremely Fast R Functions"
#>   .. .. .. .. ..$ Version           : chr "2.1.5.2"
#>   .. .. .. .. ..$ Date              : chr "2025-10-10"
#>   .. .. .. .. ..$ Authors@R         : chr "c(\n                person(\"Manos Papadakis\", email = \"rfastofficial@gmail.com\", role = c(\"aut\", \"cre\","| __truncated__
#>   .. .. .. .. ..$ Maintainer        : chr "Manos Papadakis <rfastofficial@gmail.com>"
#>   .. .. .. .. ..$ Depends           : chr "R (>= 3.5.0), Rcpp (>= 0.12.3), zigg, RcppParallel"
#>   .. .. .. .. ..$ LinkingTo         : chr "Rcpp (>= 0.12.3), RcppArmadillo, RcppParallel, zigg"
#>   .. .. .. .. ..$ Suggests          : chr "philentropy"
#>   .. .. .. .. ..$ SystemRequirements: chr "C++17"
#>   .. .. .. .. ..$ BugReports        : chr "https://github.com/RfastOfficial/Rfast/issues"
#>   .. .. .. .. ..$ URL               : chr "https://github.com/RfastOfficial/Rfast"
#>   .. .. .. .. ..$ Description       : chr "A collection of fast (utility) functions for data analysis. Column and row wise means, medians, variances, mini"| __truncated__
#>   .. .. .. .. ..$ License           : chr "GPL (>= 2.0)"
#>   .. .. .. .. ..$ NeedsCompilation  : chr "yes"
#>   .. .. .. .. ..$ Packaged          : chr "2025-10-10 12:13:26 UTC; papad"
#>   .. .. .. .. ..$ Author            : chr "Manos Papadakis [aut, cre, cph],\n  Michail Tsagris [aut],\n  Marios Dimitriadis [aut],\n  Stefanos Fafalios [a"| __truncated__
#>   .. .. .. .. ..$ Repository        : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication  : chr "2025-10-10 13:40:02 UTC"
#>   .. .. .. .. ..$ Encoding          : chr "UTF-8"
#>   .. .. .. .. ..$ Built             : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-18 05:44:39 UTC; unix"
#>   .. .. .. .. ..$ RemoteType        : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef      : chr "Rfast"
#>   .. .. .. .. ..$ RemoteRef         : chr "Rfast"
#>   .. .. .. .. ..$ RemoteRepos       : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha         : chr "2.1.5.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/Rfast/Meta/package.rds"
#>   .. .. .. ..$ spacetime   :List of 28
#>   .. .. .. .. ..$ Package          : chr "spacetime"
#>   .. .. .. .. ..$ Version          : chr "1.3-3"
#>   .. .. .. .. ..$ Title            : chr "Classes and Methods for Spatio-Temporal Data"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Edzer\",  \"Pebesma\", role = c(\"aut\", \"cre\"), email = \"edzer.pebesma@uni-muenster.de\", commen"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.0.0)"
#>   .. .. .. .. ..$ Imports          : chr "graphics, utils, stats, methods, lattice, sp (>= 1.1-0), zoo\n(>= 1.7-9), xts (>= 0.8-8), intervals"
#>   .. .. .. .. ..$ Suggests         : chr "adehabitatLT, cshapes (>= 2.0), foreign, googleVis, gstat (>=\n1.0-16), maps, mapdata, plm, raster, RColorBrewe"| __truncated__
#>   .. .. .. .. ..$ LazyData         : chr "no"
#>   .. .. .. .. ..$ Description      : chr "Classes and methods for spatio-temporal data, including space-time regular lattices, sparse lattices, irregular"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/edzer/spacetime"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/edzer/spacetime/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr"
#>   .. .. .. .. ..$ Collate          : chr "Class-xts.R Class-ST.R Class-STFDF.R Class-STSDF.R\nClass-STIDF.R Class-STTDF.R interval.R endtime.R ST-methods"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2025-02-13 15:47:01 UTC; edzer"
#>   .. .. .. .. ..$ Author           : chr "Edzer Pebesma [aut, cre] (<https://orcid.org/0000-0001-8049-7069>),\n  Benedikt Graeler [ctb],\n  Tom Gottfried"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Edzer Pebesma <edzer.pebesma@uni-muenster.de>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-02-13 16:30:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2025-04-05 03:27:15 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "spacetime"
#>   .. .. .. .. ..$ RemoteRef        : chr "spacetime"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.3-3"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/spacetime/Meta/package.rds"
#>   .. .. .. ..$ R6          :List of 27
#>   .. .. .. .. ..$ Package                : chr "R6"
#>   .. .. .. .. ..$ Title                  : chr "Encapsulated Classes with Reference Semantics"
#>   .. .. .. .. ..$ Version                : chr "2.6.1"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Winston\", \"Chang\", , \"winston@posit.co\", role = c(\"aut\", \"cre\")),\n    person(\"Posit"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Creates classes with reference semantics, similar to R's\n    built-in reference classes. Compared to reference"| __truncated__
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://r6.r-lib.org, https://github.com/r-lib/R6"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/r-lib/R6/issues"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 3.6)"
#>   .. .. .. .. ..$ Suggests               : chr "lobstr, testthat (>= 3.0.0)"
#>   .. .. .. .. ..$ Config/Needs/website   : chr "tidyverse/tidytemplate, ggplot2, microbenchmark,\nscales"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.2"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2025-02-14 21:15:19 UTC; winston"
#>   .. .. .. .. ..$ Author                 : chr "Winston Chang [aut, cre],\n  Posit Software, PBC [cph, fnd]"
#>   .. .. .. .. ..$ Maintainer             : chr "Winston Chang <winston@posit.co>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2025-02-15 00:50:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2025-04-05 02:55:50 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "R6"
#>   .. .. .. .. ..$ RemoteRef              : chr "R6"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "2.6.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/R6/Meta/package.rds"
#>   .. .. .. ..$ zoo         :List of 25
#>   .. .. .. .. ..$ Package          : chr "zoo"
#>   .. .. .. .. ..$ Version          : chr "1.8-15"
#>   .. .. .. .. ..$ Date             : chr "2025-12-15"
#>   .. .. .. .. ..$ Title            : chr "S3 Infrastructure for Regular and Irregular Time Series (Z's\nOrdered Observations)"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(given = \"Achim\", family = \"Zeileis\", role = c(\"aut\", \"cre\"), email = \"Achim.Zeileis@R-project"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "An S3 class with methods for totally ordered indexed\n             observations. It is particularly aimed at ir"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.1.0), stats"
#>   .. .. .. .. ..$ Suggests         : chr "AER, coda, chron, ggplot2 (>= 3.5.0), mondate, scales,\nstinepack, strucchange, timeDate, timeSeries, tinyplot,"| __truncated__
#>   .. .. .. .. ..$ Imports          : chr "utils, graphics, grDevices, lattice (>= 0.20-27)"
#>   .. .. .. .. ..$ License          : chr "GPL-2 | GPL-3"
#>   .. .. .. .. ..$ URL              : chr "https://zoo.R-Forge.R-project.org/"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2025-12-15 12:01:09 UTC; zeileis"
#>   .. .. .. .. ..$ Author           : chr "Achim Zeileis [aut, cre] (ORCID:\n    <https://orcid.org/0000-0003-0918-3766>),\n  Gabor Grothendieck [aut],\n "| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Achim Zeileis <Achim.Zeileis@R-project.org>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-12-15 14:30:02 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-12-16 05:16:08 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "zoo"
#>   .. .. .. .. ..$ RemoteRef        : chr "zoo"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.8-15"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/zoo/Meta/package.rds"
#>   .. .. .. ..$ lifecycle   :List of 29
#>   .. .. .. .. ..$ Package                : chr "lifecycle"
#>   .. .. .. .. ..$ Title                  : chr "Manage the Life Cycle of your Package Functions"
#>   .. .. .. .. ..$ Version                : chr "1.0.5"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Lionel\", \"Henry\", , \"lionel@posit.co\", role = c(\"aut\", \"cre\")),\n    person(\"Hadley\"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Manage the life cycle of your exported functions with shared\n    conventions, documentation badges, and user-f"| __truncated__
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://lifecycle.r-lib.org/, https://github.com/r-lib/lifecycle"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/r-lib/lifecycle/issues"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 3.6)"
#>   .. .. .. .. ..$ Imports                : chr "cli (>= 3.4.0), rlang (>= 1.1.0)"
#>   .. .. .. .. ..$ Suggests               : chr "covr, knitr, lintr (>= 3.1.0), rmarkdown, testthat (>=\n3.0.1), tibble, tidyverse, tools, vctrs, withr, xml2"
#>   .. .. .. .. ..$ VignetteBuilder        : chr "knitr"
#>   .. .. .. .. ..$ Config/Needs/website   : chr "tidyverse/tidytemplate, usethis"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2026-01-07 15:00:49 UTC; lionel"
#>   .. .. .. .. ..$ Author                 : chr "Lionel Henry [aut, cre],\n  Hadley Wickham [aut] (ORCID: <https://orcid.org/0000-0003-4757-117X>),\n  Posit Sof"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Lionel Henry <lionel@posit.co>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2026-01-08 08:00:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2026-01-09 04:36:12 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "lifecycle"
#>   .. .. .. .. ..$ RemoteRef              : chr "lifecycle"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "1.0.5"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/lifecycle/Meta/package.rds"
#>   .. .. .. ..$ fs          :List of 33
#>   .. .. .. .. ..$ Package                   : chr "fs"
#>   .. .. .. .. ..$ Title                     : chr "Cross-Platform File System Operations Based on 'libuv'"
#>   .. .. .. .. ..$ Version                   : chr "2.0.1"
#>   .. .. .. .. ..$ Authors@R                 : chr "c(\n    person(\"Jim\", \"Hester\", role = \"aut\"),\n    person(\"Hadley\", \"Wickham\", role = \"aut\"),\n   "| __truncated__
#>   .. .. .. .. ..$ Description               : chr "A cross-platform interface to file system operations, built\n    on top of the 'libuv' C library."
#>   .. .. .. .. ..$ License                   : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                       : chr "https://fs.r-lib.org, https://github.com/r-lib/fs"
#>   .. .. .. .. ..$ BugReports                : chr "https://github.com/r-lib/fs/issues"
#>   .. .. .. .. ..$ Depends                   : chr "R (>= 4.1)"
#>   .. .. .. .. ..$ Imports                   : chr "methods"
#>   .. .. .. .. ..$ Suggests                  : chr "covr, crayon, knitr, pillar (>= 1.0.0), rmarkdown, spelling,\ntestthat (>= 3.0.0), tibble (>= 1.1.0), vctrs (>= 0.3.0), withr"
#>   .. .. .. .. ..$ VignetteBuilder           : chr "knitr"
#>   .. .. .. .. ..$ SystemRequirements        : chr "libuv: libuv-devel (rpm) or libuv1-dev (deb).\nAlternatively to build the vendored libuv 'cmake' is required.\nGNU make."
#>   .. .. .. .. ..$ Config/Needs/website      : chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Config/testthat/edition   : chr "3"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep: chr "2025-04-23"
#>   .. .. .. .. ..$ Copyright                 : chr "file COPYRIGHTS"
#>   .. .. .. .. ..$ Encoding                  : chr "UTF-8"
#>   .. .. .. .. ..$ Language                  : chr "en-US"
#>   .. .. .. .. ..$ RoxygenNote               : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation          : chr "yes"
#>   .. .. .. .. ..$ Packaged                  : chr "2026-03-23 14:30:54 UTC; jeroen"
#>   .. .. .. .. ..$ Author                    : chr "Jim Hester [aut],\n  Hadley Wickham [aut],\n  Gábor Csárdi [aut],\n  Jeroen Ooms [cre],\n  libuv project contri"| __truncated__
#>   .. .. .. .. ..$ Maintainer                : chr "Jeroen Ooms <jeroenooms@gmail.com>"
#>   .. .. .. .. ..$ Repository                : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication          : chr "2026-03-24 05:20:18 UTC"
#>   .. .. .. .. ..$ Built                     : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-25 04:55:45 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef              : chr "fs"
#>   .. .. .. .. ..$ RemoteRef                 : chr "fs"
#>   .. .. .. .. ..$ RemoteRepos               : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform         : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                 : chr "2.0.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/fs/Meta/package.rds"
#>   .. .. .. ..$ ragg        :List of 32
#>   .. .. .. .. ..$ Type                             : chr "Package"
#>   .. .. .. .. ..$ Package                          : chr "ragg"
#>   .. .. .. .. ..$ Title                            : chr "Graphic Devices Based on AGG"
#>   .. .. .. .. ..$ Version                          : chr "1.5.2"
#>   .. .. .. .. ..$ Authors@R                        : chr "c(\n    person(\"Thomas Lin\", \"Pedersen\", , \"thomas.pedersen@posit.co\", role = c(\"cre\", \"aut\"),\n     "| __truncated__
#>   .. .. .. .. ..$ Maintainer                       : chr "Thomas Lin Pedersen <thomas.pedersen@posit.co>"
#>   .. .. .. .. ..$ Description                      : chr "Anti-Grain Geometry (AGG) is a high-quality and\n    high-performance 2D drawing library. The 'ragg' package pr"| __truncated__
#>   .. .. .. .. ..$ License                          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                              : chr "https://ragg.r-lib.org, https://github.com/r-lib/ragg"
#>   .. .. .. .. ..$ BugReports                       : chr "https://github.com/r-lib/ragg/issues"
#>   .. .. .. .. ..$ Imports                          : chr "systemfonts (>= 1.0.3), textshaping (>= 0.3.0)"
#>   .. .. .. .. ..$ Suggests                         : chr "covr, graphics, grid, testthat (>= 3.0.0)"
#>   .. .. .. .. ..$ LinkingTo                        : chr "systemfonts, textshaping"
#>   .. .. .. .. ..$ Config/build/compilation-database: chr "true"
#>   .. .. .. .. ..$ Config/Needs/website             : chr "ggplot2, devoid, magick, bench, tidyr, ggridges,\nhexbin, sessioninfo, pkgdown, tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Config/testthat/edition          : chr "3"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep       : chr "2025-04-25"
#>   .. .. .. .. ..$ Encoding                         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                      : chr "7.3.3"
#>   .. .. .. .. ..$ SystemRequirements               : chr "freetype2, libpng, libtiff, libjpeg, libwebp,\nlibwebpmux"
#>   .. .. .. .. ..$ NeedsCompilation                 : chr "yes"
#>   .. .. .. .. ..$ Packaged                         : chr "2026-03-16 07:50:17 UTC; thomas"
#>   .. .. .. .. ..$ Author                           : chr "Thomas Lin Pedersen [cre, aut] (ORCID:\n    <https://orcid.org/0000-0002-5147-4711>),\n  Maxim Shemanarev [aut,"| __truncated__
#>   .. .. .. .. ..$ Repository                       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication                 : chr "2026-03-23 08:50:08 UTC"
#>   .. .. .. .. ..$ Built                            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-05 04:52:07 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                     : chr "ragg"
#>   .. .. .. .. ..$ RemoteRef                        : chr "ragg"
#>   .. .. .. .. ..$ RemoteRepos                      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform                : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                        : chr "1.5.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/ragg/Meta/package.rds"
#>   .. .. .. ..$ pkgconfig   :List of 24
#>   .. .. .. .. ..$ Package          : chr "pkgconfig"
#>   .. .. .. .. ..$ Title            : chr "Private Configuration for 'R' Packages"
#>   .. .. .. .. ..$ Version          : chr "2.0.3"
#>   .. .. .. .. ..$ Author           : chr "Gábor Csárdi"
#>   .. .. .. .. ..$ Maintainer       : chr "Gábor Csárdi <csardi.gabor@gmail.com>"
#>   .. .. .. .. ..$ Description      : chr "Set configuration options on a per-package basis.\n    Options set by a given package only apply to that packag"| __truncated__
#>   .. .. .. .. ..$ License          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ LazyData         : chr "true"
#>   .. .. .. .. ..$ Imports          : chr "utils"
#>   .. .. .. .. ..$ Suggests         : chr "covr, testthat, disposables (>= 1.0.3)"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/r-lib/pkgconfig#readme"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/r-lib/pkgconfig/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2019-09-22 08:42:40 UTC; gaborcsardi"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2019-09-22 09:20:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2025-04-05 03:05:57 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "pkgconfig"
#>   .. .. .. .. ..$ RemoteRef        : chr "pkgconfig"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "2.0.3"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/pkgconfig/Meta/package.rds"
#>   .. .. .. ..$ desc        :List of 30
#>   .. .. .. .. ..$ Package                : chr "desc"
#>   .. .. .. .. ..$ Title                  : chr "Manipulate DESCRIPTION Files"
#>   .. .. .. .. ..$ Version                : chr "1.4.3"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"Gábor\", \"Csárdi\", , \"csardi.gabor@gmail.com\", role = c(\"aut\", \"cre\")),\n    person(\""| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Gábor Csárdi <csardi.gabor@gmail.com>"
#>   .. .. .. .. ..$ Description            : chr "Tools to read, write, create, and manipulate DESCRIPTION\n    files.  It is intended for packages that create o"| __truncated__
#>   .. .. .. .. ..$ License                : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                    : chr "https://desc.r-lib.org/, https://github.com/r-lib/desc"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/r-lib/desc/issues"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 3.4)"
#>   .. .. .. .. ..$ Imports                : chr "cli, R6, utils"
#>   .. .. .. .. ..$ Suggests               : chr "callr, covr, gh, spelling, testthat, whoami, withr"
#>   .. .. .. .. ..$ Config/Needs/website   : chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ Language               : chr "en-US"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.2.3"
#>   .. .. .. .. ..$ Collate                : chr "'assertions.R' 'authors-at-r.R' 'built.R' 'classes.R'\n'collate.R' 'constants.R' 'deps.R' 'desc-package.R'\n'de"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2023-12-10 11:07:50 UTC; gaborcsardi"
#>   .. .. .. .. ..$ Author                 : chr "Gábor Csárdi [aut, cre],\n  Kirill Müller [aut],\n  Jim Hester [aut],\n  Maëlle Salmon [ctb] (<https://orcid.or"| __truncated__
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2023-12-10 11:40:08 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2025-04-05 02:57:24 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "desc"
#>   .. .. .. .. ..$ RemoteRef              : chr "desc"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "1.4.3"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/desc/Meta/package.rds"
#>   .. .. .. ..$ progressr   :List of 28
#>   .. .. .. .. ..$ Package          : chr "progressr"
#>   .. .. .. .. ..$ Version          : chr "0.19.0"
#>   .. .. .. .. ..$ Title            : chr "An Inclusive, Unifying API for Progress Updates"
#>   .. .. .. .. ..$ Description      : chr "A minimal, unifying API for scripts and packages to report progress updates from anywhere including when using "| __truncated__
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Henrik\", \"Bengtsson\",\n                    role = c(\"aut\", \"cre\", \"cph\"),\n                "| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 3)"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.5.0)"
#>   .. .. .. .. ..$ Imports          : chr "digest, utils"
#>   .. .. .. .. ..$ Suggests         : chr "graphics, tcltk, beepr, cli, crayon, pbmcapply, progress,\npurrr, foreach, plyr, doFuture, future, future.apply"| __truncated__
#>   .. .. .. .. ..$ VignetteBuilder  : chr "progressr"
#>   .. .. .. .. ..$ Language         : chr "en-US"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ URL              : chr "https://progressr.futureverse.org,\nhttps://github.com/futureverse/progressr"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/futureverse/progressr/issues"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-31 06:41:07 UTC; hb"
#>   .. .. .. .. ..$ Author           : chr "Henrik Bengtsson [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0002-7579-5165>)"
#>   .. .. .. .. ..$ Maintainer       : chr "Henrik Bengtsson <henrikb@braju.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-31 08:10:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2026-04-01 04:41:36 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "progressr"
#>   .. .. .. .. ..$ RemoteRef        : chr "progressr"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.19.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/progressr/Meta/package.rds"
#>   .. .. .. ..$ pkgdown     :List of 34
#>   .. .. .. .. ..$ Package                    : chr "pkgdown"
#>   .. .. .. .. ..$ Title                      : chr "Make Static HTML Documentation for a Package"
#>   .. .. .. .. ..$ Version                    : chr "2.2.0"
#>   .. .. .. .. ..$ Authors@R                  : chr "c(\n    person(\"Hadley\", \"Wickham\", , \"hadley@posit.co\", role = c(\"aut\", \"cre\"),\n           comment "| __truncated__
#>   .. .. .. .. ..$ Description                : chr "Generate an attractive and useful website from a source\n    package.  'pkgdown' converts your documentation, v"| __truncated__
#>   .. .. .. .. ..$ License                    : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                        : chr "https://pkgdown.r-lib.org/, https://github.com/r-lib/pkgdown"
#>   .. .. .. .. ..$ BugReports                 : chr "https://github.com/r-lib/pkgdown/issues"
#>   .. .. .. .. ..$ Depends                    : chr "R (>= 4.1)"
#>   .. .. .. .. ..$ Imports                    : chr "bslib (>= 0.5.1), callr (>= 3.7.3), cli (>= 3.6.1), desc (>=\n1.4.0), downlit (>= 0.4.4), fontawesome, fs (>= 1"| __truncated__
#>   .. .. .. .. ..$ Suggests                   : chr "covr, diffviewer, evaluate (>= 0.24.0), gert, gt, htmltools,\nhtmlwidgets, knitr (>= 1.50), magick, methods, pk"| __truncated__
#>   .. .. .. .. ..$ VignetteBuilder            : chr "knitr, quarto"
#>   .. .. .. .. ..$ Config/Needs/website       : chr "usethis, servr"
#>   .. .. .. .. ..$ Config/potools/style       : chr "explicit"
#>   .. .. .. .. ..$ Config/testthat/edition    : chr "3"
#>   .. .. .. .. ..$ Config/testthat/parallel   : chr "true"
#>   .. .. .. .. ..$ Config/testthat/start-first: chr "build-article, build-quarto-article,\nbuild-reference, build"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep : chr "2025-09-07"
#>   .. .. .. .. ..$ Encoding                   : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                : chr "7.3.3"
#>   .. .. .. .. ..$ SystemRequirements         : chr "pandoc"
#>   .. .. .. .. ..$ NeedsCompilation           : chr "no"
#>   .. .. .. .. ..$ Packaged                   : chr "2025-11-05 16:48:47 UTC; hadleywickham"
#>   .. .. .. .. ..$ Author                     : chr "Hadley Wickham [aut, cre] (ORCID:\n    <https://orcid.org/0000-0003-4757-117X>),\n  Jay Hesselberth [aut] (ORCI"| __truncated__
#>   .. .. .. .. ..$ Maintainer                 : chr "Hadley Wickham <hadley@posit.co>"
#>   .. .. .. .. ..$ Repository                 : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication           : chr "2025-11-06 06:10:20 UTC"
#>   .. .. .. .. ..$ Built                      : chr "R 4.5.0; ; 2025-11-07 04:27:59 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                 : chr "any"
#>   .. .. .. .. ..$ RemotePkgRef               : chr "any::pkgdown"
#>   .. .. .. .. ..$ RemoteRef                  : chr "any::pkgdown"
#>   .. .. .. .. ..$ RemoteRepos                : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform          : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                  : chr "2.2.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/pkgdown/Meta/package.rds"
#>   .. .. .. ..$ RcppParallel:List of 28
#>   .. .. .. .. ..$ Package           : chr "RcppParallel"
#>   .. .. .. .. ..$ Type              : chr "Package"
#>   .. .. .. .. ..$ Title             : chr "Parallel Programming Tools for 'Rcpp'"
#>   .. .. .. .. ..$ Version           : chr "5.1.11-2"
#>   .. .. .. .. ..$ Authors@R         : chr "c(\n    person(\"JJ\", \"Allaire\", role = c(\"aut\"), email = \"jj@rstudio.com\"),\n    person(\"Romain\", \"F"| __truncated__
#>   .. .. .. .. ..$ Description       : chr "High level functions for parallel programming with 'Rcpp'.\n    For example, the 'parallelFor()' function can b"| __truncated__
#>   .. .. .. .. ..$ Depends           : chr "R (>= 3.0.2)"
#>   .. .. .. .. ..$ Suggests          : chr "Rcpp, RUnit, knitr, rmarkdown"
#>   .. .. .. .. ..$ SystemRequirements: chr "GNU make, Intel TBB, Windows: cmd.exe and\ncscript.exe, Solaris: g++ is required"
#>   .. .. .. .. ..$ License           : chr "GPL (>= 3)"
#>   .. .. .. .. ..$ URL               : chr "https://rcppcore.github.io/RcppParallel/,\nhttps://github.com/RcppCore/RcppParallel"
#>   .. .. .. .. ..$ BugReports        : chr "https://github.com/RcppCore/RcppParallel/issues"
#>   .. .. .. .. ..$ Biarch            : chr "TRUE"
#>   .. .. .. .. ..$ RoxygenNote       : chr "7.1.1"
#>   .. .. .. .. ..$ Encoding          : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation  : chr "yes"
#>   .. .. .. .. ..$ Packaged          : chr "2026-03-03 18:36:15 UTC; kevin"
#>   .. .. .. .. ..$ Maintainer        : chr "Kevin Ushey <kevin@rstudio.com>"
#>   .. .. .. .. ..$ Author            : chr "JJ Allaire [aut],\n  Romain Francois [aut, cph],\n  Kevin Ushey [aut, cre],\n  Gregory Vandenbrouck [aut],\n  M"| __truncated__
#>   .. .. .. .. ..$ Repository        : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication  : chr "2026-03-05 12:20:09 UTC"
#>   .. .. .. .. ..$ Built             : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-06 05:40:03 UTC; unix"
#>   .. .. .. .. ..$ RemoteType        : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef      : chr "RcppParallel"
#>   .. .. .. .. ..$ RemoteRef         : chr "RcppParallel"
#>   .. .. .. .. ..$ RemoteRepos       : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha         : chr "5.1.11-2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/RcppParallel/Meta/package.rds"
#>   .. .. .. ..$ bslib       :List of 33
#>   .. .. .. .. ..$ Package                    : chr "bslib"
#>   .. .. .. .. ..$ Title                      : chr "Custom 'Bootstrap' 'Sass' Themes for 'shiny' and 'rmarkdown'"
#>   .. .. .. .. ..$ Version                    : chr "0.10.0"
#>   .. .. .. .. ..$ Authors@R                  : chr "c(\n    person(\"Carson\", \"Sievert\", , \"carson@posit.co\", role = c(\"aut\", \"cre\"),\n           comment "| __truncated__
#>   .. .. .. .. ..$ Description                : chr "Simplifies custom 'CSS' styling of both 'shiny' and\n    'rmarkdown' via 'Bootstrap' 'Sass'. Supports 'Bootstra"| __truncated__
#>   .. .. .. .. ..$ License                    : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                        : chr "https://rstudio.github.io/bslib/, https://github.com/rstudio/bslib"
#>   .. .. .. .. ..$ BugReports                 : chr "https://github.com/rstudio/bslib/issues"
#>   .. .. .. .. ..$ Depends                    : chr "R (>= 2.10)"
#>   .. .. .. .. ..$ Imports                    : chr "base64enc, cachem, fastmap (>= 1.1.1), grDevices, htmltools\n(>= 0.5.8), jquerylib (>= 0.1.3), jsonlite, lifecy"| __truncated__
#>   .. .. .. .. ..$ Suggests                   : chr "brand.yml, bsicons, curl, fontawesome, future, ggplot2,\nknitr, lattice, magrittr, rappdirs, rmarkdown (>= 2.7)"| __truncated__
#>   .. .. .. .. ..$ Config/Needs/deploy        : chr "BH, chiflights22, colourpicker, commonmark, cpp11,\ncpsievert/chiflights22, cpsievert/histoslider, dplyr, DT,\n"| __truncated__
#>   .. .. .. .. ..$ Config/Needs/routine       : chr "chromote, desc, renv"
#>   .. .. .. .. ..$ Config/Needs/website       : chr "brio, crosstalk, dplyr, DT, ggplot2, glue,\nhtmlwidgets, leaflet, lorem, palmerpenguins, plotly, purrr,\nrprojr"| __truncated__
#>   .. .. .. .. ..$ Config/testthat/edition    : chr "3"
#>   .. .. .. .. ..$ Config/testthat/parallel   : chr "true"
#>   .. .. .. .. ..$ Config/testthat/start-first: chr "zzzz-bs-sass, fonts, zzz-precompile,\ntheme-*, rmd-*"
#>   .. .. .. .. ..$ Encoding                   : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                : chr "7.3.3"
#>   .. .. .. .. ..$ Collate                    : chr "'accordion.R' 'breakpoints.R' 'bs-current-theme.R'\n'bs-dependencies.R' 'bs-global.R' 'bs-remove.R'\n'bs-theme-"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation           : chr "no"
#>   .. .. .. .. ..$ Packaged                   : chr "2026-01-22 00:32:06 UTC; garrick"
#>   .. .. .. .. ..$ Author                     : chr "Carson Sievert [aut, cre] (ORCID:\n    <https://orcid.org/0000-0002-4958-2844>),\n  Joe Cheng [aut],\n  Garrick"| __truncated__
#>   .. .. .. .. ..$ Maintainer                 : chr "Carson Sievert <carson@posit.co>"
#>   .. .. .. .. ..$ Repository                 : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication           : chr "2026-01-26 08:10:02 UTC"
#>   .. .. .. .. ..$ Built                      : chr "R 4.5.0; ; 2026-01-27 04:33:53 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                 : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef               : chr "bslib"
#>   .. .. .. .. ..$ RemoteRef                  : chr "bslib"
#>   .. .. .. .. ..$ RemoteRepos                : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform          : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                  : chr "0.10.0"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/bslib/Meta/package.rds"
#>   .. .. .. ..$ data.table  :List of 27
#>   .. .. .. .. ..$ Package          : chr "data.table"
#>   .. .. .. .. ..$ Version          : chr "1.18.2.1"
#>   .. .. .. .. ..$ Title            : chr "Extension of `data.frame`"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.4.0)"
#>   .. .. .. .. ..$ Imports          : chr "methods"
#>   .. .. .. .. ..$ Suggests         : chr "bit64 (>= 4.0.0), bit (>= 4.0.4), R.utils (>= 2.13.0), xts,\nzoo (>= 1.8-1), yaml, knitr, markdown"
#>   .. .. .. .. ..$ Description      : chr "Fast aggregation of large data (e.g. 100GB in RAM), fast ordered joins, fast add/modify/delete of columns by gr"| __truncated__
#>   .. .. .. .. ..$ License          : chr "MPL-2.0 | file LICENSE"
#>   .. .. .. .. ..$ URL              : chr "https://r-datatable.com, https://Rdatatable.gitlab.io/data.table,\nhttps://github.com/Rdatatable/data.table"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/Rdatatable/data.table/issues"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ ByteCompile      : chr "TRUE"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n  person(\"Tyson\",\"Barrett\",        role=c(\"aut\",\"cre\"), email=\"t.barrett88@gmail.com\", comment = "| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-01-25 04:12:25 UTC; tysonbarrett"
#>   .. .. .. .. ..$ Author           : chr "Tyson Barrett [aut, cre] (ORCID:\n    <https://orcid.org/0000-0002-2137-1391>),\n  Matt Dowle [aut],\n  Arun Sr"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Tyson Barrett <t.barrett88@gmail.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-01-27 11:40:15 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-01-28 04:44:25 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "data.table"
#>   .. .. .. .. ..$ RemoteRef        : chr "data.table"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.18.2.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/data.table/Meta/package.rds"
#>   .. .. .. ..$ Rcpp        :List of 29
#>   .. .. .. .. ..$ Package          : chr "Rcpp"
#>   .. .. .. .. ..$ Title            : chr "Seamless R and C++ Integration"
#>   .. .. .. .. ..$ Version          : chr "1.1.1"
#>   .. .. .. .. ..$ Date             : chr "2026-01-07"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Dirk\", \"Eddelbuettel\", role = c(\"aut\", \"cre\"), email = \"edd@debian.org\",\n                 "| __truncated__
#>   .. .. .. .. ..$ Description      : chr "The 'Rcpp' package provides R functions as well as C++ classes which\n offer a seamless integration of R and C+"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.5.0)"
#>   .. .. .. .. ..$ Imports          : chr "methods, utils"
#>   .. .. .. .. ..$ Suggests         : chr "tinytest, inline, rbenchmark, pkgKitten (>= 0.1.2)"
#>   .. .. .. .. ..$ URL              : chr "https://www.rcpp.org,\nhttps://dirk.eddelbuettel.com/code/rcpp.html,\nhttps://github.com/RcppCore/Rcpp"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/RcppCore/Rcpp/issues"
#>   .. .. .. .. ..$ MailingList      : chr "rcpp-devel@lists.r-forge.r-project.org"
#>   .. .. .. .. ..$ RoxygenNote      : chr "6.1.1"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "Rcpp"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-01-08 14:33:45 UTC; edd"
#>   .. .. .. .. ..$ Author           : chr "Dirk Eddelbuettel [aut, cre] (ORCID:\n    <https://orcid.org/0000-0001-6419-907X>),\n  Romain Francois [aut] (O"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Dirk Eddelbuettel <edd@debian.org>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-01-10 09:50:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-01-11 07:33:45 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "Rcpp"
#>   .. .. .. .. ..$ RemoteRef        : chr "Rcpp"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.1.1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/Rcpp/Meta/package.rds"
#>   .. .. .. ..$ gstat       :List of 26
#>   .. .. .. .. ..$ Package          : chr "gstat"
#>   .. .. .. .. ..$ Version          : chr "2.1-6"
#>   .. .. .. .. ..$ Title            : chr "Spatial and Spatio-Temporal Geostatistical Modelling, Prediction\nand Simulation"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(given = \"Edzer\", \n        family = \"Pebesma\", \n        role = c(\"aut\", \"cre\"),\n        emai"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Variogram modelling; simple, ordinary and universal point or block (co)kriging; spatio-temporal kriging; sequen"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 2.10)"
#>   .. .. .. .. ..$ Imports          : chr "utils, stats, graphics, methods, lattice, sp (>= 0.9-72), zoo,\nsf (>= 0.7-2), sftime, spacetime (>= 1.2-8), stars, FNN"
#>   .. .. .. .. ..$ Suggests         : chr "fields, maps, mapdata, xts, raster, future, future.apply,\nRColorBrewer, geoR, ggplot2"
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2.0)"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/r-spatial/gstat/,\nhttps://r-spatial.github.io/gstat/"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/r-spatial/gstat/issues/"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ RoxygenNote      : chr "6.1.1"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-29 20:06:42 UTC; edzer"
#>   .. .. .. .. ..$ Author           : chr "Edzer Pebesma [aut, cre] (ORCID:\n    <https://orcid.org/0000-0001-8049-7069>),\n  Benedikt Graeler [aut]"
#>   .. .. .. .. ..$ Maintainer       : chr "Edzer Pebesma <edzer.pebesma@uni-muenster.de>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-30 05:40:07 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-30 17:14:28 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "gstat"
#>   .. .. .. .. ..$ RemoteRef        : chr "gstat"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "2.1-6"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/gstat/Meta/package.rds"
#>   .. .. .. ..$ systemfonts :List of 33
#>   .. .. .. .. ..$ Type                             : chr "Package"
#>   .. .. .. .. ..$ Package                          : chr "systemfonts"
#>   .. .. .. .. ..$ Title                            : chr "System Native Font Finding"
#>   .. .. .. .. ..$ Version                          : chr "1.3.2"
#>   .. .. .. .. ..$ Authors@R                        : chr "c(\n    person(\"Thomas Lin\", \"Pedersen\", , \"thomas.pedersen@posit.co\", role = c(\"aut\", \"cre\"),\n     "| __truncated__
#>   .. .. .. .. ..$ Description                      : chr "Provides system native access to the font catalogue. As font\n    handling varies between systems it is difficu"| __truncated__
#>   .. .. .. .. ..$ License                          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL                              : chr "https://github.com/r-lib/systemfonts,\nhttps://systemfonts.r-lib.org"
#>   .. .. .. .. ..$ BugReports                       : chr "https://github.com/r-lib/systemfonts/issues"
#>   .. .. .. .. ..$ Depends                          : chr "R (>= 3.2.0)"
#>   .. .. .. .. ..$ Imports                          : chr "base64enc, grid, jsonlite, lifecycle, tools, utils"
#>   .. .. .. .. ..$ Suggests                         : chr "covr, farver, ggplot2, graphics, knitr, ragg, rmarkdown,\nsvglite, testthat (>= 2.1.0)"
#>   .. .. .. .. ..$ LinkingTo                        : chr "cpp11 (>= 0.2.1)"
#>   .. .. .. .. ..$ VignetteBuilder                  : chr "knitr"
#>   .. .. .. .. ..$ Config/build/compilation-database: chr "true"
#>   .. .. .. .. ..$ Config/Needs/website             : chr "tidyverse/tidytemplate"
#>   .. .. .. .. ..$ Config/usethis/last-upkeep       : chr "2025-04-23"
#>   .. .. .. .. ..$ Encoding                         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                      : chr "7.3.2"
#>   .. .. .. .. ..$ SystemRequirements               : chr "fontconfig, freetype2"
#>   .. .. .. .. ..$ NeedsCompilation                 : chr "yes"
#>   .. .. .. .. ..$ Packaged                         : chr "2026-03-05 12:24:52 UTC; thomas"
#>   .. .. .. .. ..$ Author                           : chr "Thomas Lin Pedersen [aut, cre] (ORCID:\n    <https://orcid.org/0000-0002-5147-4711>),\n  Jeroen Ooms [aut] (ORC"| __truncated__
#>   .. .. .. .. ..$ Maintainer                       : chr "Thomas Lin Pedersen <thomas.pedersen@posit.co>"
#>   .. .. .. .. ..$ Repository                       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication                 : chr "2026-03-05 13:10:02 UTC"
#>   .. .. .. .. ..$ Built                            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-05 04:47:16 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                     : chr "systemfonts"
#>   .. .. .. .. ..$ RemoteRef                        : chr "systemfonts"
#>   .. .. .. .. ..$ RemoteRepos                      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform                : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                        : chr "1.3.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/systemfonts/Meta/package.rds"
#>   .. .. .. ..$ DEoptimR    :List of 27
#>   .. .. .. .. ..$ Package                         : chr "DEoptimR"
#>   .. .. .. .. ..$ Version                         : chr "1.1-4"
#>   .. .. .. .. ..$ Date                            : chr "2025-07-27"
#>   .. .. .. .. ..$ Title                           : chr "Differential Evolution Optimization in Pure R"
#>   .. .. .. .. ..$ Authors@R                       : chr "c(\n    person(c(\"Eduardo\", \"L. T.\"), \"Conceicao\",\n           role = c(\"aut\", \"cre\"),\n           em"| __truncated__
#>   .. .. .. .. ..$ URL                             : chr "svn://svn.r-forge.r-project.org/svnroot/robustbase/pkg/DEoptimR"
#>   .. .. .. .. ..$ Description                     : chr "Differential Evolution (DE) stochastic heuristic algorithms for\n    global optimization of problems with and w"| __truncated__
#>   .. .. .. .. ..$ Imports                         : chr "stats"
#>   .. .. .. .. ..$ Enhances                        : chr "robustbase"
#>   .. .. .. .. ..$ License                         : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ Author                          : chr "Eduardo L. T. Conceicao [aut, cre],\n  Martin Maechler [ctb] (ORCID: <https://orcid.org/0000-0002-8685-9910>)"
#>   .. .. .. .. ..$ Maintainer                      : chr "Eduardo L. T. Conceicao <mail@eduardoconceicao.org>"
#>   .. .. .. .. ..$ Repository                      : chr "RSPM"
#>   .. .. .. .. ..$ Repository/R-Forge/Project      : chr "robustbase"
#>   .. .. .. .. ..$ Repository/R-Forge/Revision     : chr "1009"
#>   .. .. .. .. ..$ Repository/R-Forge/DateTimeStamp: chr "2025-07-27 18:16:52"
#>   .. .. .. .. ..$ Date/Publication                : chr "2025-07-27 19:10:02 UTC"
#>   .. .. .. .. ..$ NeedsCompilation                : chr "no"
#>   .. .. .. .. ..$ Packaged                        : chr "2025-07-27 18:25:06 UTC; rforge"
#>   .. .. .. .. ..$ Encoding                        : chr "UTF-8"
#>   .. .. .. .. ..$ Built                           : chr "R 4.5.0; ; 2025-07-28 10:41:10 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                      : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                    : chr "DEoptimR"
#>   .. .. .. .. ..$ RemoteRef                       : chr "DEoptimR"
#>   .. .. .. .. ..$ RemoteRepos                     : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform               : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                       : chr "1.1-4"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/DEoptimR/Meta/package.rds"
#>   .. .. .. ..$ xfun        :List of 28
#>   .. .. .. .. ..$ Package          : chr "xfun"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "Supporting Functions for Packages Maintained by 'Yihui Xie'"
#>   .. .. .. .. ..$ Version          : chr "0.57"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n  person(\"Yihui\", \"Xie\", role = c(\"aut\", \"cre\", \"cph\"), email = \"xie@yihui.name\", comment = c(O"| __truncated__
#>   .. .. .. .. ..$ Description      : chr "Miscellaneous functions commonly used in other packages maintained by 'Yihui Xie'."
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.2.0)"
#>   .. .. .. .. ..$ Imports          : chr "grDevices, stats, tools"
#>   .. .. .. .. ..$ Suggests         : chr "testit, parallel, codetools, methods, rstudioapi, tinytex (>=\n0.30), mime, litedown (>= 0.6), commonmark, knit"| __truncated__
#>   .. .. .. .. ..$ License          : chr "MIT + file LICENSE"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/yihui/xfun"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/yihui/xfun/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote      : chr "7.3.3"
#>   .. .. .. .. ..$ VignetteBuilder  : chr "litedown"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-03-19 17:27:02 UTC; yihui"
#>   .. .. .. .. ..$ Author           : chr "Yihui Xie [aut, cre, cph] (ORCID:\n    <https://orcid.org/0000-0003-0645-5666>, URL: https://yihui.org),\n  Wus"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Yihui Xie <xie@yihui.name>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-03-20 16:50:02 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-21 04:43:41 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "xfun"
#>   .. .. .. .. ..$ RemoteRef        : chr "xfun"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.57"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/xfun/Meta/package.rds"
#>   .. .. .. ..$ knitr       :List of 30
#>   .. .. .. .. ..$ Package           : chr "knitr"
#>   .. .. .. .. ..$ Type              : chr "Package"
#>   .. .. .. .. ..$ Title             : chr "A General-Purpose Package for Dynamic Report Generation in R"
#>   .. .. .. .. ..$ Version           : chr "1.51"
#>   .. .. .. .. ..$ Authors@R         : chr "c(\n    person(\"Yihui\", \"Xie\", role = c(\"aut\", \"cre\"), email = \"xie@yihui.name\", comment = c(ORCID = "| __truncated__
#>   .. .. .. .. ..$ Description       : chr "Provides a general-purpose tool for dynamic report generation in R\n    using Literate Programming techniques."
#>   .. .. .. .. ..$ Depends           : chr "R (>= 3.6.0)"
#>   .. .. .. .. ..$ Imports           : chr "evaluate (>= 0.15), highr (>= 0.11), methods, tools, xfun (>=\n0.52), yaml (>= 2.1.19)"
#>   .. .. .. .. ..$ Suggests          : chr "bslib, DBI (>= 0.4-1), digest, formatR, gifski, gridSVG,\nhtmlwidgets (>= 0.7), jpeg, JuliaCall (>= 0.11.1), ma"| __truncated__
#>   .. .. .. .. ..$ License           : chr "GPL"
#>   .. .. .. .. ..$ URL               : chr "https://yihui.org/knitr/"
#>   .. .. .. .. ..$ BugReports        : chr "https://github.com/yihui/knitr/issues"
#>   .. .. .. .. ..$ Encoding          : chr "UTF-8"
#>   .. .. .. .. ..$ VignetteBuilder   : chr "litedown, knitr"
#>   .. .. .. .. ..$ SystemRequirements: chr "Package vignettes based on R Markdown v2 or\nreStructuredText require Pandoc (http://pandoc.org). The\nfunction"| __truncated__
#>   .. .. .. .. ..$ Collate           : chr "'block.R' 'cache.R' 'citation.R' 'hooks-html.R' 'plot.R'\n'utils.R' 'defaults.R' 'concordance.R' 'engine.R' 'hi"| __truncated__
#>   .. .. .. .. ..$ RoxygenNote       : chr "7.3.3"
#>   .. .. .. .. ..$ NeedsCompilation  : chr "no"
#>   .. .. .. .. ..$ Packaged          : chr "2025-12-19 15:19:25 UTC; yihui"
#>   .. .. .. .. ..$ Author            : chr "Yihui Xie [aut, cre] (ORCID: <https://orcid.org/0000-0003-0645-5666>,\n    URL: https://yihui.org),\n  Abhranee"| __truncated__
#>   .. .. .. .. ..$ Maintainer        : chr "Yihui Xie <xie@yihui.name>"
#>   .. .. .. .. ..$ Repository        : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication  : chr "2025-12-20 14:30:02 UTC"
#>   .. .. .. .. ..$ Built             : chr "R 4.5.0; ; 2025-12-21 04:20:33 UTC; unix"
#>   .. .. .. .. ..$ RemoteType        : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef      : chr "knitr"
#>   .. .. .. .. ..$ RemoteRef         : chr "knitr"
#>   .. .. .. .. ..$ RemoteRepos       : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha         : chr "1.51"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/knitr/Meta/package.rds"
#>   .. .. .. ..$ nlme        :List of 24
#>   .. .. .. .. ..$ Package         : chr "nlme"
#>   .. .. .. .. ..$ Version         : chr "3.1-168"
#>   .. .. .. .. ..$ Date            : chr "2025-03-31"
#>   .. .. .. .. ..$ Priority        : chr "recommended"
#>   .. .. .. .. ..$ Title           : chr "Linear and Nonlinear Mixed Effects Models"
#>   .. .. .. .. ..$ Authors@R       : chr "c(person(\"José\", \"Pinheiro\", role = \"aut\", comment = \"S version\"),\n             person(\"Douglas\", \""| __truncated__
#>   .. .. .. .. ..$ Contact         : chr "see 'MailingList'"
#>   .. .. .. .. ..$ Description     : chr "Fit and compare Gaussian linear and nonlinear mixed-effects models."
#>   .. .. .. .. ..$ Depends         : chr "R (>= 3.6.0)"
#>   .. .. .. .. ..$ Imports         : chr "graphics, stats, utils, lattice"
#>   .. .. .. .. ..$ Suggests        : chr "MASS, SASmixed"
#>   .. .. .. .. ..$ LazyData        : chr "yes"
#>   .. .. .. .. ..$ Encoding        : chr "UTF-8"
#>   .. .. .. .. ..$ License         : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ BugReports      : chr "https://bugs.r-project.org"
#>   .. .. .. .. ..$ MailingList     : chr "R-help@r-project.org"
#>   .. .. .. .. ..$ URL             : chr "https://svn.r-project.org/R-packages/trunk/nlme/"
#>   .. .. .. .. ..$ NeedsCompilation: chr "yes"
#>   .. .. .. .. ..$ Packaged        : chr "2025-03-31 11:19:09 UTC; ripley"
#>   .. .. .. .. ..$ Author          : chr "José Pinheiro [aut] (S version),\n  Douglas Bates [aut] (up to 2007),\n  Saikat DebRoy [ctb] (up to 2002),\n  D"| __truncated__
#>   .. .. .. .. ..$ Maintainer      : chr "R Core Team <R-core@R-project.org>"
#>   .. .. .. .. ..$ Repository      : chr "CRAN"
#>   .. .. .. .. ..$ Date/Publication: chr "2025-03-31 11:21:01 UTC"
#>   .. .. .. .. ..$ Built           : chr "R 4.5.3; x86_64-pc-linux-gnu; 2026-03-11 09:35:33 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/nlme/Meta/package.rds"
#>   .. .. .. ..$ htmltools   :List of 31
#>   .. .. .. .. ..$ Type                : chr "Package"
#>   .. .. .. .. ..$ Package             : chr "htmltools"
#>   .. .. .. .. ..$ Title               : chr "Tools for HTML"
#>   .. .. .. .. ..$ Version             : chr "0.5.9"
#>   .. .. .. .. ..$ Authors@R           : chr "c(\n    person(\"Joe\", \"Cheng\", , \"joe@posit.co\", role = \"aut\"),\n    person(\"Carson\", \"Sievert\", , "| __truncated__
#>   .. .. .. .. ..$ Description         : chr "Tools for HTML generation and output."
#>   .. .. .. .. ..$ License             : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL                 : chr "https://github.com/rstudio/htmltools,\nhttps://rstudio.github.io/htmltools/"
#>   .. .. .. .. ..$ BugReports          : chr "https://github.com/rstudio/htmltools/issues"
#>   .. .. .. .. ..$ Depends             : chr "R (>= 2.14.1)"
#>   .. .. .. .. ..$ Imports             : chr "base64enc, digest, fastmap (>= 1.1.0), grDevices, rlang (>=\n1.0.0), utils"
#>   .. .. .. .. ..$ Suggests            : chr "Cairo, markdown, ragg, shiny, testthat, withr"
#>   .. .. .. .. ..$ Enhances            : chr "knitr"
#>   .. .. .. .. ..$ Config/Needs/check  : chr "knitr"
#>   .. .. .. .. ..$ Config/Needs/website: chr "rstudio/quillt, bench"
#>   .. .. .. .. ..$ Encoding            : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote         : chr "7.3.3"
#>   .. .. .. .. ..$ Collate             : chr "'colors.R' 'fill.R' 'html_dependency.R' 'html_escape.R'\n'html_print.R' 'htmltools-package.R' 'images.R' 'known"| __truncated__
#>   .. .. .. .. ..$ NeedsCompilation    : chr "yes"
#>   .. .. .. .. ..$ Packaged            : chr "2025-12-02 23:32:55 UTC; cpsievert"
#>   .. .. .. .. ..$ Author              : chr "Joe Cheng [aut],\n  Carson Sievert [aut, cre] (ORCID:\n    <https://orcid.org/0000-0002-4958-2844>),\n  Barret "| __truncated__
#>   .. .. .. .. ..$ Maintainer          : chr "Carson Sievert <carson@posit.co>"
#>   .. .. .. .. ..$ Repository          : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication    : chr "2025-12-04 11:50:09 UTC"
#>   .. .. .. .. ..$ Built               : chr "R 4.5.0; x86_64-pc-linux-gnu; 2025-12-04 16:32:21 UTC; unix"
#>   .. .. .. .. ..$ RemoteType          : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef        : chr "htmltools"
#>   .. .. .. .. ..$ RemoteRef           : chr "htmltools"
#>   .. .. .. .. ..$ RemoteRepos         : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform   : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha           : chr "0.5.9"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/htmltools/Meta/package.rds"
#>   .. .. .. ..$ igraph      :List of 39
#>   .. .. .. .. ..$ Package                            : chr "igraph"
#>   .. .. .. .. ..$ Title                              : chr "Network Analysis and Visualization"
#>   .. .. .. .. ..$ Version                            : chr "2.2.3"
#>   .. .. .. .. ..$ Authors@R                          : chr "c(\n    person(\"Gábor\", \"Csárdi\", , \"csardi.gabor@gmail.com\", role = \"aut\",\n           comment = c(ORC"| __truncated__
#>   .. .. .. .. ..$ Description                        : chr "Routines for simple graphs and network analysis. It can\n    handle large graphs very well and provides functio"| __truncated__
#>   .. .. .. .. ..$ License                            : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL                                : chr "https://r.igraph.org/, https://igraph.org/,\nhttps://igraph.discourse.group/"
#>   .. .. .. .. ..$ BugReports                         : chr "https://github.com/igraph/rigraph/issues"
#>   .. .. .. .. ..$ Depends                            : chr "methods, R (>= 3.5.0)"
#>   .. .. .. .. ..$ Imports                            : chr "cli, graphics, grDevices, lifecycle, magrittr, Matrix,\npkgconfig (>= 2.0.0), rlang (>= 1.1.0), stats, utils, vctrs"
#>   .. .. .. .. ..$ Suggests                           : chr "ape (>= 5.7-0.1), callr, decor, digest, igraphdata, knitr,\nrgl (>= 1.3.14), rmarkdown, scales, stats4, tcltk, "| __truncated__
#>   .. .. .. .. ..$ Enhances                           : chr "graph"
#>   .. .. .. .. ..$ LinkingTo                          : chr "cpp11 (>= 0.5.0)"
#>   .. .. .. .. ..$ VignetteBuilder                    : chr "knitr"
#>   .. .. .. .. ..$ Config/build/compilation-database  : chr "false"
#>   .. .. .. .. ..$ Config/build/never-clean           : chr "true"
#>   .. .. .. .. ..$ Config/comment/compilation-database: chr "Generate manually with\npkgload:::generate_db() for faster pkgload::load_all()"
#>   .. .. .. .. ..$ Config/Needs/build                 : chr "r-lib/roxygen2, devtools, irlba, pkgconfig,\nigraph/igraph.r2cdocs, moodymudskipper/devtag"
#>   .. .. .. .. ..$ Config/Needs/coverage              : chr "covr"
#>   .. .. .. .. ..$ Config/Needs/website               : chr "here, readr, tibble, xmlparsedata, xml2"
#>   .. .. .. .. ..$ Config/testthat/edition            : chr "3"
#>   .. .. .. .. ..$ Config/testthat/parallel           : chr "true"
#>   .. .. .. .. ..$ Config/testthat/start-first        : chr "aaa-auto, vs-es, scan, vs-operators,\nweakref, watts.strogatz.game"
#>   .. .. .. .. ..$ Encoding                           : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote                        : chr "7.3.3.9000"
#>   .. .. .. .. ..$ SystemRequirements                 : chr "libxml2 (optional), glpk (>= 4.57, optional)"
#>   .. .. .. .. ..$ NeedsCompilation                   : chr "yes"
#>   .. .. .. .. ..$ Packaged                           : chr "2026-04-07 03:24:11 UTC; kirill"
#>   .. .. .. .. ..$ Author                             : chr "Gábor Csárdi [aut] (ORCID: <https://orcid.org/0000-0001-7098-9676>),\n  Tamás Nepusz [aut] (ORCID: <https://orc"| __truncated__
#>   .. .. .. .. ..$ Maintainer                         : chr "Kirill Müller <kirill@cynkra.com>"
#>   .. .. .. .. ..$ Repository                         : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication                   : chr "2026-04-07 06:40:02 UTC"
#>   .. .. .. .. ..$ Built                              : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-04-08 04:45:01 UTC; unix"
#>   .. .. .. .. ..$ RemoteType                         : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef                       : chr "igraph"
#>   .. .. .. .. ..$ RemoteRef                          : chr "igraph"
#>   .. .. .. .. ..$ RemoteRepos                        : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform                  : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha                          : chr "2.2.3"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/igraph/Meta/package.rds"
#>   .. .. .. ..$ rmarkdown   :List of 31
#>   .. .. .. .. ..$ Type                   : chr "Package"
#>   .. .. .. .. ..$ Package                : chr "rmarkdown"
#>   .. .. .. .. ..$ Title                  : chr "Dynamic Documents for R"
#>   .. .. .. .. ..$ Version                : chr "2.31"
#>   .. .. .. .. ..$ Authors@R              : chr "c(\n    person(\"JJ\", \"Allaire\", , \"jj@posit.co\", role = \"aut\"),\n    person(\"Yihui\", \"Xie\", , \"xie"| __truncated__
#>   .. .. .. .. ..$ Description            : chr "Convert R Markdown documents into a variety of formats."
#>   .. .. .. .. ..$ License                : chr "GPL-3"
#>   .. .. .. .. ..$ URL                    : chr "https://github.com/rstudio/rmarkdown,\nhttps://pkgs.rstudio.com/rmarkdown/"
#>   .. .. .. .. ..$ BugReports             : chr "https://github.com/rstudio/rmarkdown/issues"
#>   .. .. .. .. ..$ Depends                : chr "R (>= 3.0)"
#>   .. .. .. .. ..$ Imports                : chr "bslib (>= 0.2.5.1), evaluate (>= 0.13), fontawesome (>=\n0.5.0), htmltools (>= 0.5.1), jquerylib, jsonlite, kni"| __truncated__
#>   .. .. .. .. ..$ Suggests               : chr "digest, dygraphs, fs, rsconnect, downlit (>= 0.4.0), katex\n(>= 1.4.0), sass (>= 0.4.0), shiny (>= 1.6.0), test"| __truncated__
#>   .. .. .. .. ..$ VignetteBuilder        : chr "knitr"
#>   .. .. .. .. ..$ Config/Needs/website   : chr "rstudio/quillt, pkgdown"
#>   .. .. .. .. ..$ Config/testthat/edition: chr "3"
#>   .. .. .. .. ..$ Encoding               : chr "UTF-8"
#>   .. .. .. .. ..$ RoxygenNote            : chr "7.3.3"
#>   .. .. .. .. ..$ SystemRequirements     : chr "pandoc (>= 1.14) - http://pandoc.org"
#>   .. .. .. .. ..$ NeedsCompilation       : chr "no"
#>   .. .. .. .. ..$ Packaged               : chr "2026-03-25 18:23:38 UTC; yihui"
#>   .. .. .. .. ..$ Author                 : chr "JJ Allaire [aut],\n  Yihui Xie [aut, cre] (ORCID: <https://orcid.org/0000-0003-0645-5666>),\n  Christophe Dervi"| __truncated__
#>   .. .. .. .. ..$ Maintainer             : chr "Yihui Xie <xie@yihui.name>"
#>   .. .. .. .. ..$ Repository             : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication       : chr "2026-03-26 16:30:02 UTC"
#>   .. .. .. .. ..$ Built                  : chr "R 4.5.0; ; 2026-03-27 05:26:54 UTC; unix"
#>   .. .. .. .. ..$ RemoteType             : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef           : chr "rmarkdown"
#>   .. .. .. .. ..$ RemoteRef              : chr "rmarkdown"
#>   .. .. .. .. ..$ RemoteRepos            : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform      : chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha              : chr "2.31"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/rmarkdown/Meta/package.rds"
#>   .. .. .. ..$ xts         :List of 28
#>   .. .. .. .. ..$ Package          : chr "xts"
#>   .. .. .. .. ..$ Type             : chr "Package"
#>   .. .. .. .. ..$ Title            : chr "eXtensible Time Series"
#>   .. .. .. .. ..$ Version          : chr "0.14.2"
#>   .. .. .. .. ..$ Authors@R        : chr "c(\n  person(given=c(\"Jeffrey\",\"A.\"), family=\"Ryan\", role=c(\"aut\",\"cph\")),\n  person(given=c(\"Joshua"| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.6.0), zoo (>= 1.7-12)"
#>   .. .. .. .. ..$ Imports          : chr "methods"
#>   .. .. .. .. ..$ LinkingTo        : chr "zoo"
#>   .. .. .. .. ..$ Suggests         : chr "timeSeries, timeDate, tseries, chron, tinytest"
#>   .. .. .. .. ..$ LazyLoad         : chr "yes"
#>   .. .. .. .. ..$ Description      : chr "Provide for uniform handling of R's different time-based data classes by extending zoo, maximizing native forma"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL              : chr "https://joshuaulrich.github.io/xts/,\nhttps://github.com/joshuaulrich/xts"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/joshuaulrich/xts/issues"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-02-27 23:01:11 UTC; josh"
#>   .. .. .. .. ..$ Author           : chr "Jeffrey A. Ryan [aut, cph],\n  Joshua M. Ulrich [cre, aut],\n  Ross Bennett [ctb],\n  Corwin Joy [ctb]"
#>   .. .. .. .. ..$ Maintainer       : chr "Joshua M. Ulrich <josh.m.ulrich@gmail.com>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-02-28 06:10:03 UTC"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-03-01 04:37:38 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "xts"
#>   .. .. .. .. ..$ RemoteRef        : chr "xts"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "0.14.2"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/xts/Meta/package.rds"
#>   .. .. .. ..$ ismev       :List of 22
#>   .. .. .. .. ..$ Package          : chr "ismev"
#>   .. .. .. .. ..$ Version          : chr "1.43"
#>   .. .. .. .. ..$ Date             : chr "2025-08-12"
#>   .. .. .. .. ..$ Title            : chr "An Introduction to Statistical Modeling of Extreme Values"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Eric\", \"Gilleland\", role = \"cre\", \n                  email = \"eric.gilleland@colostate.edu\","| __truncated__
#>   .. .. .. .. ..$ Author           : chr "Eric Gilleland [cre] (ORCID: <https://orcid.org/0000-0002-8058-7643>),\n  Janet E. Heffernan [aut],\n  Alec G. "| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Eric Gilleland <eric.gilleland@colostate.edu>"
#>   .. .. .. .. ..$ Depends          : chr "R (>= 2.10.0), mgcv"
#>   .. .. .. .. ..$ Description      : chr "Functions to support the computations carried out in\n  `An Introduction to Statistical Modeling of Extreme Val"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ NeedsCompilation : chr "no"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Packaged         : chr "2025-08-12 13:52:09 UTC; gille"
#>   .. .. .. .. ..$ Date/Publication : chr "2025-08-20 17:00:15 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; ; 2025-08-21 10:42:30 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "ismev"
#>   .. .. .. .. ..$ RemoteRef        : chr "ismev"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "1.43"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/ismev/Meta/package.rds"
#>   .. .. .. ..$ compiler    :List of 10
#>   .. .. .. .. ..$ Package    : chr "compiler"
#>   .. .. .. .. ..$ Version    : chr "4.5.3"
#>   .. .. .. .. ..$ Priority   : chr "base"
#>   .. .. .. .. ..$ Title      : chr "The R Compiler Package"
#>   .. .. .. .. ..$ Author     : chr "Luke Tierney <luke-tierney@uiowa.edu>"
#>   .. .. .. .. ..$ Maintainer : chr "R Core Team <do-use-Contact-address@r-project.org>"
#>   .. .. .. .. ..$ Contact    : chr "R-help mailing list <r-help@r-project.org>"
#>   .. .. .. .. ..$ Description: chr "Byte code compiler for R."
#>   .. .. .. .. ..$ License    : chr "Part of R 4.5.3"
#>   .. .. .. .. ..$ Built      : chr "R 4.5.3; ; 2026-03-11 09:30:23 UTC; unix"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/opt/R/4.5.3/lib/R/library/compiler/Meta/package.rds"
#>   .. .. .. ..$ sp          :List of 27
#>   .. .. .. .. ..$ Package          : chr "sp"
#>   .. .. .. .. ..$ Version          : chr "2.2-1"
#>   .. .. .. .. ..$ Title            : chr "Classes and Methods for Spatial Data"
#>   .. .. .. .. ..$ Authors@R        : chr "c(person(\"Edzer\", \"Pebesma\", role = c(\"aut\", \"cre\"),\n\t\t\t\temail = \"edzer.pebesma@uni-muenster.de\""| __truncated__
#>   .. .. .. .. ..$ Depends          : chr "R (>= 3.5.0), methods"
#>   .. .. .. .. ..$ Imports          : chr "utils, stats, graphics, grDevices, lattice, grid"
#>   .. .. .. .. ..$ Suggests         : chr "RColorBrewer, gstat, deldir, knitr, maps, mapview, rmarkdown,\nsf, terra, raster"
#>   .. .. .. .. ..$ Description      : chr "Classes and methods for spatial\n  data; the classes document where the spatial location information\n  resides"| __truncated__
#>   .. .. .. .. ..$ License          : chr "GPL (>= 2)"
#>   .. .. .. .. ..$ URL              : chr "https://github.com/edzer/sp/ https://edzer.github.io/sp/"
#>   .. .. .. .. ..$ BugReports       : chr "https://github.com/edzer/sp/issues"
#>   .. .. .. .. ..$ Collate          : chr "bpy.colors.R AAA.R Class-CRS.R CRS-methods.R Class-Spatial.R\nSpatial-methods.R projected.R Class-SpatialPoints"| __truncated__
#>   .. .. .. .. ..$ VignetteBuilder  : chr "knitr"
#>   .. .. .. .. ..$ NeedsCompilation : chr "yes"
#>   .. .. .. .. ..$ Packaged         : chr "2026-02-12 22:13:59 UTC; edzer"
#>   .. .. .. .. ..$ Author           : chr "Edzer Pebesma [aut, cre],\n  Roger Bivand [aut],\n  Barry Rowlingson [ctb],\n  Virgilio Gomez-Rubio [ctb],\n  R"| __truncated__
#>   .. .. .. .. ..$ Maintainer       : chr "Edzer Pebesma <edzer.pebesma@uni-muenster.de>"
#>   .. .. .. .. ..$ Repository       : chr "RSPM"
#>   .. .. .. .. ..$ Date/Publication : chr "2026-02-13 10:30:02 UTC"
#>   .. .. .. .. ..$ Encoding         : chr "UTF-8"
#>   .. .. .. .. ..$ Built            : chr "R 4.5.0; x86_64-pc-linux-gnu; 2026-02-15 04:37:28 UTC; unix"
#>   .. .. .. .. ..$ RemoteType       : chr "standard"
#>   .. .. .. .. ..$ RemotePkgRef     : chr "sp"
#>   .. .. .. .. ..$ RemoteRef        : chr "sp"
#>   .. .. .. .. ..$ RemoteRepos      : chr "https://packagemanager.posit.co/cran/__linux__/noble/latest"
#>   .. .. .. .. ..$ RemotePkgPlatform: chr "x86_64-pc-linux-gnu-ubuntu-24.04"
#>   .. .. .. .. ..$ RemoteSha        : chr "2.2-1"
#>   .. .. .. .. ..- attr(*, "class")= chr "packageDescription"
#>   .. .. .. .. ..- attr(*, "file")= chr "/home/runner/work/_temp/Library/sp/Meta/package.rds"
#>   .. .. ..$ matprod    : chr "default"
#>   .. .. ..$ BLAS       : chr "/usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3"
#>   .. .. ..$ LAPACK     : chr "/usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so"
#>   .. .. ..$ LA_version : chr "3.12.0"
#>   .. .. ..- attr(*, "class")= chr "sessionInfo"
#>   .. ..$ run.time: 'difftime' num 0.380061864852905
#>   .. .. ..- attr(*, "units")= chr "secs"
#>  - attr(*, "plR.track")= chr "000"
#>  - attr(*, "class")= chr "plR"
```

## Input File Formats

polylinkR requires three input files:

### 1. ObjInfo.txt (Gene Information)

Contains gene-level information:

- **objID**: Unique gene identifier
- **objStat**: Gene-level statistic (e.g., p-value, z-score)
- Optional columns: chr, startpos, endpos (for coordinate-based
  analysis)

### 2. SetInfo.txt (Gene Set Information)

Contains gene set definitions:

- **setID**: Unique set identifier
- **setName**: Human-readable set name
- **setSource**: Source database (optional)

### 3. SetObj.txt (Gene-Set Mappings)

Contains gene-to-set relationships:

- **setID**: Set identifier (must match SetInfo.txt)
- **objID**: Gene identifier (must match ObjInfo.txt)

## Important Notes

### Active Development

This package is under active development. Results should be
independently validated before use in production or clinical decisions.

### Algorithm Constraints

The core statistical algorithms (permutation nulls, GPD fitting,
rescaling, pruning) are frozen and should not be modified without
explicit maintainer approval.

## Next Steps

After reading the data, you would typically proceed through the
pipeline:

1.  **plR_permute()** - Generate permutation-based null distributions
2.  **plR_rescale()** - Account for genetic autocorrelation
3.  **plR_prune()** - Identify significant gene sets

See the function documentation
([`?plR_permute`](../reference/plR_permute.md),
[`?plR_rescale`](../reference/plR_rescale.md),
[`?plR_prune`](../reference/plR_prune.md)) for details on subsequent
steps.

## Session Information

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] polylinkR_0.5.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] dqrng_0.4.1           sass_0.4.10           future_1.70.0        
#>  [4] robustbase_0.99-7     distances_0.1.13      lattice_0.22-9       
#>  [7] listenv_0.10.1        digest_0.6.39         magrittr_2.0.5       
#> [10] evaluate_1.0.5        grid_4.5.3            iterators_1.0.14     
#> [13] fastmap_1.2.0         Matrix_1.7-4          foreach_1.5.2        
#> [16] jsonlite_2.0.0        mgcv_1.9-4            codetools_0.2-20     
#> [19] textshaping_1.0.5     jquerylib_0.1.4       cli_3.6.6            
#> [22] rlang_1.2.0           zigg_0.0.2            parallelly_1.46.1    
#> [25] future.apply_1.20.2   splines_4.5.3         intervals_0.15.5     
#> [28] cachem_1.1.0          yaml_2.3.12           FNN_1.1.4.1          
#> [31] tools_4.5.3           parallel_4.5.3        doFuture_1.2.1       
#> [34] globals_0.19.1        Rfast_2.1.5.2         spacetime_1.3-3      
#> [37] R6_2.6.1              zoo_1.8-15            lifecycle_1.0.5      
#> [40] fs_2.0.1              ragg_1.5.2            pkgconfig_2.0.3      
#> [43] desc_1.4.3            progressr_0.19.0      pkgdown_2.2.0        
#> [46] RcppParallel_5.1.11-2 bslib_0.10.0          data.table_1.18.2.1  
#> [49] Rcpp_1.1.1            gstat_2.1-6           systemfonts_1.3.2    
#> [52] DEoptimR_1.1-4        xfun_0.57             knitr_1.51           
#> [55] nlme_3.1-168          htmltools_0.5.9       igraph_2.2.3         
#> [58] rmarkdown_2.31        xts_0.14.2            ismev_1.43           
#> [61] compiler_4.5.3        sp_2.2-1
```
