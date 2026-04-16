# Print method for `plR` objects

Prints an overview of a `plR` object, including its processing history
and available contents.

## Usage

``` r
# S3 method for class 'plR'
print(x, ...)
```

## Arguments

- x:

  `plR` class object; typically output from a `polylinkR` function.

- ...:

  Additional arguments passed to `print`.

## Value

An overview of the `plR` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Provide a snapshot of core files (assuming `my_plR` is a valid plR object)
print(my_plR)

# Simply typing the object name also prints the snapshot
my_plR
} # }
```
