# Batch LEiDA pipeline across subjects

Apply
[`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md)
to multiple subjects in one call. Accepts either a named list of \[N x
Tmax\] matrices or a 3-D array \[N x Tmax x subjects\].

## Usage

``` r
batch_leida(data, ...)
```

## Arguments

- data:

  Named list of \[N x Tmax\] numeric matrices, **or** a 3-D numeric
  array with dimensions \[N, Tmax, subjects\]. For arrays, subject names
  are taken from `dimnames(data)[[3]]` if set, otherwise `"sub1"`,
  `"sub2"`, etc.

- ...:

  Additional arguments passed to
  [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md)
  (e.g. `flp`, `fhi`, `delt`, `filter`).

## Value

Named list of `dynR_leida` objects, one per subject.

## See also

[`stack_leida()`](https://dynr.circadia-lab.uk/reference/stack_leida.md),
[`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md)

## Examples

``` r
set.seed(1)
ts_list <- list(
  sub01 = matrix(rnorm(10 * 200), nrow = 10),
  sub02 = matrix(rnorm(10 * 200), nrow = 10)
)
results <- batch_leida(ts_list, filter = FALSE)
names(results)
#> [1] "sub01" "sub02"
```
