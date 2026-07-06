# Plot Kuramoto synchrony time series

Line plot of the Kuramoto order parameter R(t) over time, with a dashed
horizontal line at the mean. Uses the dynR periwinkle/brick-red palette.

## Usage

``` r
plot_synchrony(synchrony, title = "Kuramoto order parameter")
```

## Arguments

- synchrony:

  Numeric vector. Kuramoto R(t), as returned by
  [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md) or
  the `$synchrony` element of
  [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md).

- title:

  Character. Plot title.

## Value

A `ggplot` object.

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10)
ph  <- hilbert_phases(ts)
kop <- kuramoto(ph)
plot_synchrony(kop$synchrony)
```
