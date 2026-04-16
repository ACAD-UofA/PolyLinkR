# Gene set enrichment on scores rescaled for genetic autocorrelation

Rescales gene set scores to account for dependencies within gene sets
based on autocorrelated genes and re-evaluates gene set enrichment on
the rescaled scores. Note that the `obj.info` file must contain
appropriately labelled genomic coordinates for autocorrelation to be
estimated and gene set rescaling to be performed, otherwise the function
will exit with an error message.

## Usage

``` r
plR_rescale(
  plR.input,
  rescale = TRUE,
  fast = TRUE,
  ac = NULL,
  cgm.bin = 30,
  cgm.range = "auto",
  cgm.wt.max = 0.05,
  emp.bayes = "auto",
  min.rho = 1e-05,
  verbose = TRUE,
  n.cores = "auto",
  fut.plan = "auto"
)
```

## Arguments

- plR.input:

  `plR` class object; output from
  [`polylinkR::plR_permute`](plR_permute.md). Required.

- rescale:

  `logical`; should gene set score (i.e.,`setScore.std`) rescaling be
  performed? Defaults to `TRUE`. If `FALSE`, only inter-gene
  autocorrelation is estimated without revising enrichment testing for
  rescaled (i.e., decorrelated) gene set scores. This is useful for
  exploring how parameter settings impact gene set autocorrelation. The
  rescaling step can be performed later by passing the output to
  `plR_rescale`, or by providing the `ac` object to the `ac` argument.

- fast:

  `logical`; should the gene set score rescaling be performed using the
  fast mode? Defaults to `TRUE`, whereby the analytical scaling factor
  \\1 / \sqrt{1 + (m - 1)\hat{\rho}}\\ is calculated from the covariance
  matrix `ac` and used to adjust the empirical null and GPD (scale and
  threshold) parameters. If `FALSE`, the full empirical rescaling is
  performed instead, by recalculating the score for every permuted gene
  set. This is much more computationally intensive and is intended for
  users interested in comparing analytical and empirical results.

- ac:

  `data.frame` or `data.table`; user-provided object containing
  inter-gene autocorrelations. Include all unique pairs of genes with
  non-zero autocorrelation. Genes must be identified by `objID`s in
  columns `objID.A` and `objID.B`. The corresponding autocorrelation
  value must be in `gamma`. Defaults to `NULL`.

- cgm.bin:

  `numeric`; mimimum number of gene pairs required for bins of empirical
  covariances. Default value = `30`. Must be in range `[10, 1e3]`.
  Initially an exponential grid of bin sizes will be generated,
  favouring smaller bins at short distances and larger bins for more
  distant gene pairs. Bins with too few genes will be successively
  merged with the larger bins until the minimum number of genes are
  reached. This condition is evaluated across all chromosomes, ensuring
  no bin is smaller than the minimum value (other than the final bins).

- cgm.range:

  `numeric`; maximum inter-gene lag used to evaluate autocovariance.
  Defaults to `2e6` bp or `2` cM, depending on genetic distance measure.
  Must be within interval `[1e5, 5e7]` bp or `[0.1, 50]` cM.

- cgm.wt.max:

  `numeric`; maximum probability weight for a single lag window.
  Defaults to `0.1`. Set to `1` if no upper bound is desired. Note that
  the lower bound is reset to `2 / no. fitted lags` if the chosen value
  falls below this limit.

- emp.bayes:

  `numeric`; A character string indicating the empirical Bayes rescaling
  framework employed for chromosome-level covariance beta coefficients.
  Users may choose to fit either the `full` or `reduced` model, with the
  the former model also allowing for the two coefficient to be
  correlated. The default setting (`auto`) runs the `full` model if more
  than `15` chromosomes are present, otherwise the `reduced` model is
  fit. In cases when the `full` model is assessed, the `reduced` model
  is also fitted and a likelihood-based test is used to decide whether a
  correlated random-effects structure is supported.

- min.rho:

  `numeric`; minimum estimated correlation between two gene sets; values
  below this are set to `0`. Defaults to `1e-5`. Range `(0, 0.01]`.

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

Rescaled set scores (i.e., corrected for correlations between genes
within sets) are recorded in `set.info` as `setScore.rs`. Results from
the revised enrichment tests are shown in `setScore.rs.p`.

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
# Assuming `my_plr` is the result of `polylinkR::plR_permute`

# Example 1: Basic usage
new_plr <- plR_rescale(plR.input = my_plr)

# Example 2: Estimate autocorrelation only (no rescaling)
new_plr <- plR_rescale(
   plR.input = my_plr,
   rescale = FALSE
)

# Example 3: Custom variogram model and cores
new_plr <- plR_rescale(
  plR.input = my_plr,
  vg.model  = c("Exp","Sph"),
  n.cores   = 4
)

# Example 4: Using a user-provided autocovariance object
# Assuming `my_ac` is a valid data.table or data.frame
new_plr <- plR_rescale(
   plR.input = my_plr,
   ac        = my_ac
)

# Or using the `new_plr` object generated by Example 2
new_plr <- plR_rescale(plR.input = my_rescaled_plr)
} # }
```
