# Plot a functional connectivity matrix

Render a square correlation or phase-locking matrix as a heatmap using
the dynR diverging palette (deep indigo -\> periwinkle -\> brick red),
centred at zero. Matches the colour scheme used in the dynR vignettes
and pkgdown site.

## Usage

``` r
plot_fc(fc_matrix, title = "Functional connectivity", limits = c(-1, 1))
```

## Arguments

- fc_matrix:

  Numeric matrix \[N x N\]. Pearson correlations or phase-locking
  values.

- title:

  Character. Plot title. Default `"Functional connectivity"`.

- limits:

  Numeric vector `c(lo, hi)`. Colour scale limits. Default `c(-1, 1)`.

## Value

A `ggplot` object.

## Examples

``` r
data(fc, package = "dynR")
plot_fc(fc[1:20, 1:20])
```
