# Edge-centric cofluctuation analysis

Compute edge time series and root-sum-square (RSS) cofluctuations from
z-standardised BOLD signal. Implements the edge-centric FC framework of
Esfahlani et al. (2020) and Faskowitz et al. (2020).

## Usage

``` r
cofluct(timeseries, k = 1)
```

## Arguments

- timeseries:

  Numeric matrix \[N × Tmax\]. BOLD signal with N parcels as rows and
  Tmax timepoints as columns.

- k:

  Integer. Upper-triangle offset. `k = 1` (default) excludes the
  diagonal; `k = 0` includes it.

## Value

A list with:

- edge_ts:

  Numeric matrix \[n_edges × Tmax\]. Edge time series, one row per
  unique parcel pair.

- rss:

  Numeric vector \[Tmax\]. Root-sum-square cofluctuation at each
  timepoint.

## References

Esfahlani, F. Z. et al. (2020). High-amplitude cofluctuations in
cortical activity drive functional connectivity. *PNAS*, 117(45),
28393–28401.
[doi:10.1073/pnas.2005531117](https://doi.org/10.1073/pnas.2005531117)

Faskowitz, J. et al. (2020). Edge-centric functional network
representations of human cerebral cortex reveal overlapping system-level
architecture. *Nature Neuroscience*, 23(12), 1644–1654.
[doi:10.1038/s41593-020-00719-y](https://doi.org/10.1038/s41593-020-00719-y)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
res <- cofluct(ts)
dim(res$edge_ts)  # n_edges x 200
#> [1]  45 200
```
