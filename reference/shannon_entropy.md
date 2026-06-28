# Shannon entropy

Estimate Shannon entropy from a numeric vector, with optional bit-depth
discretisation (as used internally by
[`kuramoto()`](https://CoDe-Neuro.github.io/dynR/reference/kuramoto.md)).

## Usage

``` r
shannon_entropy(x, base = 2, n_bits = NULL)
```

## Arguments

- x:

  Numeric vector.

- base:

  Numeric. Logarithm base. Default is 2 (entropy in bits).

- n_bits:

  Integer or `NULL`. If supplied, values are scaled by \\2^{n\\bits}\\
  and rounded before computing entropy. Default is `NULL`.

## Value

Numeric scalar. Shannon entropy.

## Examples

``` r
x <- sample(1:4, 100, replace = TRUE)
shannon_entropy(x)
#> [1] 1.984508
```
