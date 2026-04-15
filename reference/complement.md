# Generate data with specified correlation to a vector

Generate data with specified correlation to a vector

## Usage

``` r
complement(y, corr = 0.5)
```

## Arguments

- y:

  Vector to correlate with.

- corr:

  Target correlation value (between -1 and 1).

## Value

A vector with the specified correlation to `y`.

## Examples

``` r
y <- rnorm(100)
x <- complement(y, corr = 0.8)
cor(x, y)
#> [1] 0.8
```
