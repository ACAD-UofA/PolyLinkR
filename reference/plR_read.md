# Read and validate required files for polylinkR

Reads and validates the input files for a `polylinkR` analysis. Handles
file path resolution, column name canonicalisation, data validation, and
optional filtering and data generation tasks like oordinate conversion
or gene score calculation.

## Usage

``` r
plR_read(
  input.path = NULL,
  obj.info.path = NULL,
  set.info.path = NULL,
  set.obj.path = NULL,
  var.info.path = NULL,
  rec.rate.path = NULL,
  min.set.n = 2L,
  max.set.n = Inf,
  group = NULL,
  map.fun = "kosambi",
  obj.buffer = "auto",
  obj.stat.fun = "non.param",
  bin.size = 250L,
  obj.in = NULL,
  obj.out = NULL,
  set.in = NULL,
  set.out = NULL,
  set.merge = 0.95,
  rem.genes = FALSE,
  verbose = TRUE
)
```

## Arguments

- input.path:

  `character`; path to the directory with input files. Defaults to
  `NULL`. Not required if using separate file paths. Automatically
  searches for required files with the following labels:

  - `setinfo`: a data frame with set-level information.

  - `objinfo`: a data frame with object-level information.

  - `setobj`: a data frame mapping objects to sets.

  It also searches for optional files:

  - `recrate`: an optional file for genetic coordinate conversion.

  - `varinfo`: an optional file for gene score generation.

  Allowable file names include capitalisation of each word in the label
  and use of internal separators (e.g., `set.info`, `set_info`,
  `SetInfo`, `Set.Info`, `Set_Info`). Also see the section on required
  and optional input files for the columns needed in each input file.

- obj.info.path:

  `character`; path to the `obj.info` file. Defaults to `NULL`. Ignored
  if `input.path` is provided; otherwise all required file paths must be
  specified.

- set.info.path:

  `character`; path to the `set.info` file. Defaults to `NULL`. Ignored
  if `input.path` is provided; otherwise all required file paths must be
  specified.

- set.obj.path:

  `character`; path to the `set.obj` file. Defaults to `NULL`. Ignored
  if `input.path` is provided; otherwise all required file paths must be
  specified.

- var.info.path:

  `character`; path to the `var.info` file. Defaults to `NULL`. Optional
  file used for gene score generation. Ignored if `input.path` is
  provided.

- rec.rate.path:

  `character`; path to the `rec.rate` file. Defaults to `NULL`. Optional
  file used for genetic coordinate conversion. Ignored if `input.path`
  is provided.

- min.set.n:

  `integer`; minimum size of gene sets to be retained. Defaults to `2`.
  Must be in the range `[2L, max.set.n)`.

- max.set.n:

  `integer`; maximum size of gene sets to be retained. Defaults to
  `Inf`. Must be in the range `(min.set.n, Inf)`.

- group:

  `character`; label used to identify input files within a directory.
  Defaults to `NULL`.

- map.fun:

  `character`; mapping function to convert physical to genetic
  distances. Options are `"Haldane"`, `"Kosambi"` (default),
  `"Carter-Falconer"`, and `"Morgan"`.

- obj.buffer:

  `numeric`; interval around genes (in base pairs) to include when
  assigning values from the `var.info` file. Defaults to `1e4` if
  `var.info` is provided and user does not set a value, otherwise it is
  set to 0 if score assignment is not performed. User values must be in
  the range `[0, 1e5L]`. Note that if the user provides their own gene
  scores in the `obj.info ` input file (i.e. not computed from the
  `var.info` file), then the start and end positions must include any
  buffer used to bin scores, otherwise polylinkR deconfounding and
  autocorrelation inference will not be performed appropriately.

- obj.stat.fun:

  `character`; function used to correct maximum gene scores based on the
  number of overlapping summary statistics (SNPs or windows). Default is
  `non.param`, a robust non-parametric method that uses binned data to
  calculate median and median absolute deviation (MAD) to normalise
  scores. Alternatively, `lm.logN` applies a linear regression to the
  log-transformed SNP / bin counts (assumes a roughly linear
  relationship is appropriate). In both cases, expected scores are
  estimated and gene scores calculated as the residual value. Ignored if
  no var.info file is provided.

- bin.size:

  `integer`; gene set size interval for non-parametric correction.
  Defaults to `250L`. Must be in the range `[50L, 1e3L]`; ignored if the
  parametric function is used.

- obj.in:

  `character` or `numeric` vector; `objID`s of genes to explicitly
  retain. Defaults to `NULL`.

- obj.out:

  `character` or `numeric` vector; `objID`s of genes to explicitly
  remove. Defaults to `NULL`.

- set.in:

  `character` or `numeric` vector; `setID`s of gene sets to explicitly
  retain. Defaults to `NULL`.

- set.out:

  `character` or `numeric` vector; `setID`s of gene sets to explicitly
  remove. Defaults to `NULL`.

- set.merge:

  `numeric`; minimum proportion of shared genes for merging gene sets.
  Defaults to `0.95`. Must be in the range `(0, 1]`.

- rem.genes:

  `logical`; should genes with identical genomic positions be removed?
  Defaults to `FALSE`.

- verbose:

  `logical`; should progress messages be printed to the console?
  Defaults to `TRUE`.

## Value

A `plR` S3 object containing three complementary datasets:

- `obj.info`:

  `data.table` containing information for each gene (object).

- `set.info`:

  `data.table` containing information for each gene set.

- `set.obj`:

  `data.table` containing the mapping of genes to gene sets.

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

## Required input file structure

The following three files are compulsory for a `polylinkR` analysis and
extend the format introduced in Polysel. All files must be comma-
separated (`.csv`) or tab-separated (`.tsv`), with a header (see
examples at <https://github.com/CMPG/polysel/tree/master/data>).

- `set.info`:

  A `data.table` (and `data.frame`) with gene set information.

  - **setID**: `character` or `factor` vector of unique gene set
    identifiers. Required.

- `obj.info`:

  A `data.table` (and `data.frame`) with gene (objects) information.

  - **objID**: `character` or `factor` vector of unique gene
    identifiers. Required.

  - **objStat**: optional `numeric` vector of pre-computed gene scores.
    Required if `var.info` is absent; otherwise computed from
    `var.info`.

  - **CovX**: `numeric` vector of covariate scores, where `X` is a
    positive integer denoting covariate number. Required for
    deconfounding gene scores (`objStat`) in `plR_permute`.

  - **chr**: `character` or `numeric` chromosome / contig labels.
    Required for gene score deconfounding and gene set score
    decorrelation in and `plR_permute` and `plR_rescale`, respectively.

  - **startpos**: `numeric` gene start position in base pairs. If the
    user is providing their own gene scores, then this position
    represents the original position minus any buffer used to assign
    summary statistics to genes. Required for gene score deconfounding
    and gene set score decorrelation in and `plR_permute` and
    `plR_rescale`, respectively.

  - **endpos**: `numeric` gene end position in base pairs. If the user
    is providing their own gene scores, then this position represents
    the original position plus any buffer used to assign summary
    statistics to genes. Required for gene score deconfounding and gene
    set score decorrelation in and `plR_permute` and `plR_rescale`,
    respectively.

  - **startpos.base**: `numeric` gene start position in base pairs.
    Default lower gene boundary (ignoring buffer used in gene
    assignment). Created internally if absent or gene scores are
    estimated (i.e. var.info provided). Required for gene score
    deconfounding and set score decorrelation in `plR_permute` and
    `plR_rescale`, respectively.

  - **endpos.base**: `numeric` gene end position base pairs. Default
    upper gene boundary (ignoring buffer used in gene assignment).
    Created internally if absent or gene scores are estimated (i.e.
    var.info provided). Required for gene score deconfounding and set
    score decorrelation in `plR_permute` and `plR_rescale`,
    respectively.

- `set.obj`:

  A `data.table` (and `data.frame`) mapping genes to gene sets.

  - **setID**: `character` or `factor` vector of unique gene set
    identifiers. Required.

  - **objID**: `character` or `factor` vector of unique gene
    identifiers. Required.

## Optional input file structure

These optional comma-separated (`.csv`) or tab-separated (`.tsv`) files
provide the following `plR_read` functionality:

- `var.info`:

  Contains information used to compute gene scores (`objStat`).

  - **chr**: `character` or `numeric` chromosome / contig labels.
    Required.

  - **pos**: `numeric` position where statistic was evaluated (e.g., SNP
    or central point in window). Required.

  - **value**: `numeric` statistical value from SNP / window score.
    Required.

- `rec.rate`:

  Used to transform genetic coordinates from physical to genetic
  distances. Requires `startpos` and `endpos` in `obj.info` to be base
  pair coordinates. Uses HapMap format (see labelled examples at
  <https://zenodo.org/records/11437540>).

  - **chr**: `character` or `numeric` chromosome / contig labels.
    Required.

  - **pos**: `numeric` base pair position of upstream marker for the
    recombination interval. Required.

  - **rate**: `numeric` recombination rate in cM per bp in the
    downstream interval. Required.

  - **map**: `numeric` genetic distance in centiMorgans (cM). Optional;
    will be calculated from recombination rates if absent.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Read all files from a single folder
my_plr <- plR_read(input.path = "path/to/files")

# Example 2: Read all files with "POP1" in label from a single folder
my_plr <- plR_read(
   input.path = "path/to/files",
   group      = "POP1"
)

# Example 3: Relax gene set merging criteria
my_plr <- plR_read(
   input.path = "path/to/files",
   set.merge  = 0.50
)

# Example 4: Specify separate file paths and remove user-specified sets
my_plr <- plR_read(
  set.info.path = "path/to/set.info",
  set.obj.path  = "path/to/set.obj",
  obj.info.path = "path/to/obj.info",
  set.out       = c("set1", "set2")
)

# Example 5: Generate gene scores from var.info using regression
# and include a 50 kb buffer around genes
my_plr <- plR_read(
  set.info.path = "path/to/set.info",
  set.obj.path  = "path/to/set.obj",
  obj.info.path = "path/to/obj.info",
  var.info.path = "path/to/var.info",
  obj.stat.fun   = "lm.logN", obj.buffer = 5e4
)

# Example 6: Convert distances to cM using rec.rate file
my_plr <- plR_read(
  set.info.path = "path/to/set.info",
  set.obj.path  = "path/to/set.obj",
  obj.info.path = "path/to/obj.info",
  rec.rate.path = "path/to/rec.rate"
)
} # }
```
