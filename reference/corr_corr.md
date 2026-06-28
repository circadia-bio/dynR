# Correlation of correlations matrix

Compute the correlation between edge time series across all timepoints,
producing a \[Tmax × Tmax\] "correlation of correlations" matrix (Hansen
et al., 2015).

## Usage

``` r
corr_corr(timeseries, k = 1)
```

## Arguments

- timeseries:

  Numeric matrix \[N × Tmax\]. BOLD signal.

- k:

  Integer. Upper-triangle offset passed to
  [`cofluct()`](https://dynr.circadia-lab.uk/reference/cofluct.md).
  Default 1.

## Value

Numeric matrix \[Tmax × Tmax\].

## References

Hansen, E. C. A. et al. (2015). Functional connectivity dynamics:
Modeling the switching behavior of the resting state. *NeuroImage*, 105,
525–535.
[doi:10.1016/j.neuroimage.2014.11.001](https://doi.org/10.1016/j.neuroimage.2014.11.001)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
cc <- corr_corr(ts)
dim(cc)  # 100 x 100
#> [1] 100 100
```
