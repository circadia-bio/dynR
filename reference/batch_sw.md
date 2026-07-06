# Batch sliding-window pipeline across subjects

Apply
[`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md)
to multiple subjects in one call. Accepts either a named list of \[N x
Tmax\] matrices or a 3-D array \[N x Tmax x subjects\].

## Usage

``` r
batch_sw(data, window, ...)
```

## Arguments

- data:

  Named list of \[N x Tmax\] matrices, or a 3-D array \[N, Tmax,
  subjects\]. See
  [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md).

- window:

  Integer. Window size in timepoints passed to
  [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md).

- ...:

  Additional arguments passed to
  [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md).

## Value

Named list of `dynR_sw` objects, one per subject.

## See also

[`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md),
[`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md)

## Examples

``` r
set.seed(2)
ts_list <- list(
  sub01 = matrix(rnorm(10 * 200), nrow = 10),
  sub02 = matrix(rnorm(10 * 200), nrow = 10)
)
results <- batch_sw(ts_list, window = 20, filter = FALSE)
names(results)
#> [1] "sub01" "sub02"
```
