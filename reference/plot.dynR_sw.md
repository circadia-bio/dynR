# Plot method for dynR_sw objects

Plot method for dynR_sw objects

## Usage

``` r
# S3 method for class 'dynR_sw'
plot(x, type = c("rss", "fc"), ...)
```

## Arguments

- x:

  A `dynR_sw` object as returned by
  [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md).

- type:

  Character. `"rss"` (default) plots the RSS cofluctuation time series;
  `"fc"` plots the time-averaged sliding-window FC matrix via
  [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md).

- ...:

  Additional arguments passed to the underlying plot function.

## Value

A `ggplot` object.
