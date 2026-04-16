# Summary method for `plR` objects

Provides a list of significant gene sets from the most recently applied
step of the `polylinkR` workflow.

## Usage

``` r
# S3 method for class 'plR'
summary(object, sig = 0.05, ...)
```

## Arguments

- object:

  `plR` class object; typically the output of a `polylinkR` function.

- sig:

  `numeric`; significance threshold for filtering p- or q-values. Must
  be in the range `(0, 1)`. Defaults to `0.05`.

- ...:

  Additional arguments passed to `summary`.

## Value

A summary of significant gene sets in the `plR` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage (assuming `my_plR` is a valid plR object)
summary(my_plR)

# Using a stricter significance threshold
summary(my_plR, sig = 0.01)
} # }
```
