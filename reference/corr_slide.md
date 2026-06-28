# Sliding window correlation

Compute functional connectivity matrices over sliding windows of a BOLD
timeseries (Hansen et al., 2015).

## Usage

``` r
corr_slide(timeseries, window, step = NULL)
```

## Arguments

- timeseries:

  Numeric matrix \[N × Tmax\]. BOLD signal with N parcels as rows and
  Tmax timepoints as columns.

- window:

  Integer. Window size in timepoints.

- step:

  Integer. Step between window onsets. Defaults to `window`
  (non-overlapping windows).

## Value

A list with:

- corr_mats:

  Array \[N, N, n_windows\] of Pearson correlation matrices.

- idx:

  Integer vector of 1-indexed window onset positions.

## References

Hansen, E. C. A. et al. (2015). Functional connectivity dynamics:
Modeling the switching behavior of the resting state. *NeuroImage*, 105,
525–535.
[doi:10.1016/j.neuroimage.2014.11.001](https://doi.org/10.1016/j.neuroimage.2014.11.001)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
res <- corr_slide(ts, window = 20)
dim(res$corr_mats)  # 10 x 10 x 10
#> [1] 10 10 10
```
