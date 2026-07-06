# Plot method for dynR_leida objects

Plot method for dynR_leida objects

## Usage

``` r
# S3 method for class 'dynR_leida'
plot(x, type = c("synchrony", "fc"), ...)
```

## Arguments

- x:

  A `dynR_leida` object as returned by
  [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md).

- type:

  Character. `"synchrony"` (default) plots the Kuramoto order parameter
  time series via
  [`plot_synchrony()`](https://dynr.circadia-lab.uk/reference/plot_synchrony.md);
  `"fc"` plots the time-averaged phase-locking matrix via
  [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md).

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A `ggplot` object.
