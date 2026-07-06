# Stack LEiDA eigenvectors across subjects

Combine the LEiDA eigenvector matrices from a list of `dynR_leida`
results into a single matrix (or data frame) ready for cross-subject
K-means clustering.

## Usage

``` r
stack_leida(batch_result, add_subject_id = TRUE)
```

## Arguments

- batch_result:

  Named list of `dynR_leida` objects, as returned by
  [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md).

- add_subject_id:

  Logical. Prepend a `subject` column with the subject name? Default
  `TRUE`.

## Value

If `add_subject_id = TRUE`: a data frame with columns `subject` and one
column per parcel (named `V1`, `V2`, ...). If `FALSE`: a plain numeric
matrix \[total_timepoints x N\].

## See also

[`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md)

## Examples

``` r
set.seed(3)
ts_list <- list(
  sub01 = matrix(rnorm(10 * 200), nrow = 10),
  sub02 = matrix(rnorm(10 * 200), nrow = 10)
)
results <- batch_leida(ts_list, filter = FALSE)
leida_all <- stack_leida(results)
dim(leida_all)
#> [1] 360  11
```
