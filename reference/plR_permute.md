# Gene set enrichment with control for confounding covariates

Performs gene set enrichment testing using a gene-wise permutation
procedure, while accounting for confounding variables. Note that the
deconfounding step requires that appropriately labeled covariate columns
are provided in the `obj.info` file (otherwise standard gene score
residuals are passed to permutations). Applies a Generalized Pareto
Distribution (GPD) to the tail of the empirical null for accurate
estimation of small p-values.

## Usage

``` r
plR_permute(
  plR.input,
  permute = TRUE,
  n.perm = 500000L,
  n.boot = 30L,
  alt = "upper",
  md.meth = "robust",
  kern.wt.max = 0.05,
  kern.bound = "auto",
  kern.func = "normal",
  kern.scale = "auto",
  gpd.cutoff = 5000L/n.perm,
  seed = NULL,
  verbose = TRUE,
  n.cores = "auto",
  fut.plan = "auto"
)
```

## Arguments

- plR.input:

  `plR` class object; output from [`polylinkR::plR_read`](plR_read.md).
  Required.

- permute:

  `logical`; should the permutation step be performed? Defaults to
  `TRUE`. If `FALSE`, only gene score deconfounding is performed without
  estimating set scores. If `TRUE`, users must provide confounder
  covariates in `obj.info`. This is useful for exploring how parameter
  settings impact deconfounded gene scores; the permutation step can be
  performed later by passing the output back into `plR_permute`.

- n.perm:

  `integer`; number of permutations. Defaults to `1e5`. Must be in the
  range `[5e4L, Inf)` and be exactly divisible by `1e4`.

- n.boot:

  `integer`; number of bootstrap replicates for null inference. Defaults
  to `30L`. Must be in the range `[5, Inf)`.

- alt:

  `character`; direction of the hypothesis test. Defaults to `"upper"`
  for enrichment in the upper tail (large set scores). Alternatively,
  `"lower"` tests for enrichment in the lower tail (small values). When
  `"lower"` is chosen, data are internally negated in functions
  performing p-value estimation.

- md.meth:

  `character`; determines whether raw covariate data or ranks are used
  in Mahalanobis distance calculations. Defaults to `"robust"`, where
  Mahalanobis distances are converted to ranks and Spearman's metric is
  used to calculate the covariance matrix. The alternate option `"raw"`
  uses the original covariates and Pearson's covariance for scaling.

- kern.wt.max:

  `numeric`; maximum probability weight for a single gene. Defaults to
  `0.05`. Must be in the range `(1 / (n.genes - 1), 1]`. Set to `1` if
  no upper bound is desired. Ignored if covariate columns are not
  detected in `obj.info`.

- kern.bound:

  `numeric`; flanking region around each gene where weights of
  overlapping genes inside the region are set to 0. Weights of partially
  overlapping genes are downscaled by the proportion overlapping the
  excluded region. Defaults to 0.1 Mbp or 0.1 cM (depending on genetic
  coordinates). Range `(0, Inf]`; `0` denotes standard gene boundaries
  and `Inf` excludes all genes on the same chromosome. Ignored if
  appropriate covariate columns are not detected in `obj.info`.

- kern.func:

  `character`; kernel function used to generate probability weights from
  distances between the focal gene and other genes in confounder space.
  Default is `"normal"` (Gaussian kernel). Alternate options include
  `"exponential"` and `"inverse"`. Ignored if covariate columns are not
  detected in `obj.info`.

- kern.scale:

  `numeric`; scalar used in the kernel function to convert Mahalanobis
  distances to probabilities. Defaults to `2` for Gaussian decay,
  `log(10)` for exponential decay, or `2` for inverse decay. Range
  `(0, Inf]`. Ignored if covariate columns are not detected in
  `obj.info`.

- gpd.cutoff:

  `numeric`; threshold tail probability at which to apply GPD tail
  fitting. Defaults to `500 / n.perm`. Must be in the range
  `[max(c(1e-04, 500 / n.perm)), 0.05]`. The lower bound constraint
  ensures that a minimum of `500` exceedances are available for GPD
  estimation, while also ensuring compatibility with the empirical CDF,
  where the lowest evaluated quantile is `1e-4`.

- seed:

  `integer`; random seed for reproducibility. Preserved across
  subsequent polylinkR functions (`plR_permute`, `plR_rescale`).
  Defaults to `NULL`, in which case a seed is generated automatically.
  Must be within `[-.Machine$integer.max, .Machine$integer.max]`.

- verbose:

  `logical`; should progress messages be printed to the console?
  Defaults to `TRUE`.

- n.cores:

  `integer`; number of cores for parallel processing. Defaults to `1` or
  `maximum cores - 1`. Must be in the range `[1, maximum cores]`.

- fut.plan:

  `character`; parallel backend from the `future` package. Defaults to
  user `n.cores` choice or checks available cores, choosing
  `"sequential"` on single-core and `"multisession"` on multi-core
  systems. Options: `"multisession"`, `"multicore"`, `"cluster"`,
  `"sequential"`.

## Value

A `plR` S3 object containing:

- `obj.info`:

  `data.table` for each gene (object).

- `set.info`:

  `data.table` for each gene set.

- `set.obj`:

  `data.table` mapping genes to gene sets.

Deconfounded gene scores are recorded in `obj.info` as `objStat.std`.
Set scores are recorded in `set.info` as `setScore.std`. Results from
enrichment tests are in `setScore.std.p`.

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
# Assuming `my_plr` is the result of `polylinkR::plR_read`

# Example 1: Basic usage
new_plr <- plR_permute(plR.input = my_plr)

# Example 2: Less permutations, more bootstraps
new_plr <- plR_permute(
   plR.input = my_plr,
   n.perm    = 1e5,
   n.boot    = 100
)

# Example 3: Modified covariate handling, single processor
new_plr <- plR_permute(
  plR.input = my_plr,
  kern.wt.max  = 0.2,
  md.meth   = "raw",
  n.cores   = 1
)

# Example 4: Modified GPD estimation, user-specified seed
new_plr <- plR_permute(
  plR.input  = my_plr,
  gpd.cutoff = 0.01,
  seed       = 1000
)

# Example 5: Only deconfound scores (no enrichment analysis)
new_plr <- plR_permute(
   plR.input = my_plr,
   permute =   FALSE
)

# Example 6: Use deconfounded scores from Example 5 to run enrichment
new_plr <- plR_permute(plR.input = new_plr)
} # }
```
