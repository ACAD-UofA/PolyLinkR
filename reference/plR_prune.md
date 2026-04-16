# Gene set enrichment adjusting for multiple testing and shared genes

Adjusts gene set scores to account dependencies between gene sets using
a pruning routine involving sequential p-value ranking and removal of
shared genes. Performs multiple testing correction on the pruned gene
set p-values. Note that this step require that at more than one gene set
with shared genes is present, otherwise the function will exit with an
error message.

## Usage

``` r
plR_prune(
  plR.input,
  n.fdr = 300L,
  est.pi0 = TRUE,
  tolerance = 0.001,
  verbose = TRUE,
  n.cores = "auto",
  fut.plan = "auto"
)
```

## Arguments

- plR.input:

  `plR` class object; output from
  [`polylinkR::plR_permute`](plR_permute.md) or
  [`polylinkR::plR_rescale`](plR_rescale.md). Required.

- n.fdr:

  `integer`; number of times to replicate the pruning procedure. Used to
  estimate FDR-corrected p (q) values using a histogram method. Defaults
  to `300L`. Range `[100, Inf)` and must be exactly divisible by `100`.

- est.pi0:

  `logical`; should `pi0` be estimated during FDR correction? Defaults
  to `TRUE`. If `FALSE`, `pi0` is set to `1`.

- tolerance:

  `numeric`; minimum difference between successive `pi0` estimates for
  estimation to stop. Defaults to `1e-3`. Range `(0, 0.01]`.

- verbose:

  `logical`; should progress reporting be enabled? Defaults to `TRUE`.

- n.cores:

  `integer`; number of cores for parallel processing. Defaults to `1` or
  `maximum cores - 1`. Must be in `[1, maximum cores]`.

- fut.plan:

  `character`; parallel backend from the `future` package. Defaults to
  user `n.cores` choice or checks cores, choosing `"sequential"` on
  single-core and `"multisession"` on multi-core systems. Options:
  `"multisession"`, `"multicore"`, `"cluster"`, `"sequential"`.

## Value

A `plR` S3 object containing:

- `obj.info`:

  `data.table` for each gene (object).

- `set.info`:

  `data.table` for each gene set.

- `set.obj`:

  `data.table` mapping genes to gene sets.

Pruned set scores (i.e., corrected for between–gene-set correlations
caused by shared genes) are recorded in `set.info` as `setScore.pr`.
Corrected scores are shown in `setScore.pr.p`. Results from multiple
testing correction are shown in `setScore.pr.q`.

All `plR` S3 objects include auxiliary information and summaries as
attributes:

- `plr.data`:

  Reusable datasets and parameters, including GPD fitting results and
  autocovariance estimation.

- `plr.args`:

  Argument settings used in the function.

- `plr.session`:

  R session information and function run time.

- `plr.track`:

  Internal tracking number indicating functional steps.

- `class`:

  S3 object class (i.e., `plr`).

Each attribute can be accessed using
[`attr()`](https://rdrr.io/r/base/attr.html) or
[`attributes()`](https://rdrr.io/r/base/attributes.html). The first four
attribute classes aggregate information over successive functions. For
example, to access the `plr.data` attribute for the `plR` object output
after running `plR_permute`, use `attributes(X)$plR.data$permute.data`,
where `X` is the object name. Similarly, the arguments used in
`plR_read` are in `attributes(X)$plR.args$read.args`.

The primary data structure of the `plR` object can be accessed using
[`print()`](https://rdrr.io/r/base/print.html) or by simply typing the
object's name.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming `my_plr` is the result of `polylinkR::plR_permute` or
# `polylinkR::plR_rescale`

# Example 1: Basic usage
new_plr <- plR_prune(plR.input = my_plr)

# Example 2: Increase iterations for null p-value distribution
new_plr <- plR_prune(
   plR.input = my_plr,
   n.fdr = 1000L
)

# Example 3: Do not estimate pi0 (sets pi0 = 1)
new_plr <- plR_prune(
   plR.input = my_plr,
   est.pi0 = FALSE
)
} # }
```
