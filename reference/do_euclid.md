# Euclidean distance between consecutive points

Compute the Euclidean distance between each consecutive pair of rows in
a matrix — typically a trajectory through PC space.

## Usage

``` r
do_euclid(x)
```

## Arguments

- x:

  Numeric matrix \[n_points × n_dims\].

## Value

Numeric vector of length `nrow(x)`. The first element is always 0 (no
previous point).

## Examples

``` r
set.seed(1)
pcs <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
d <- do_euclid(pcs)
length(d)  # 50
#> [1] 50
```
