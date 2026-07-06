# Stack synchrony vectors across subjects

Combine the Kuramoto synchrony time series from a list of `dynR_leida`
results into a tidy data frame with one row per timepoint per subject.

## Usage

``` r
stack_synchrony(batch_result)
```

## Arguments

- batch_result:

  Named list of `dynR_leida` objects.

## Value

A data frame with columns:

- subject:

  Subject identifier.

- timepoint:

  Timepoint index (within subject).

- synchrony:

  Kuramoto R(t).

- metastability:

  Subject-level metastability (constant within subject).

## See also

[`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md)

## Examples

``` r
set.seed(4)
ts_list <- list(
  sub01 = matrix(rnorm(10 * 200), nrow = 10),
  sub02 = matrix(rnorm(10 * 200), nrow = 10)
)
results <- batch_leida(ts_list, filter = FALSE)
df_sync <- stack_synchrony(results)
head(df_sync)
#>         subject timepoint  synchrony metastability
#> sub01.1   sub01         1 0.18511952     0.1481606
#> sub01.2   sub01         2 0.14466531     0.1481606
#> sub01.3   sub01         3 0.05097755     0.1481606
#> sub01.4   sub01         4 0.18215059     0.1481606
#> sub01.5   sub01         5 0.18066710     0.1481606
#> sub01.6   sub01         6 0.13160931     0.1481606
```
