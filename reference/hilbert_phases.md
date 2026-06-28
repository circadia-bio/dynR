# Hilbert transform phase extraction

Compute the instantaneous phase time series for all parcels/voxels via
the analytic signal (Hilbert transform). Each parcel timeseries is
demeaned before transformation.

## Usage

``` r
hilbert_phases(timeseries)
```

## Arguments

- timeseries:

  Numeric matrix \[N × Tmax\]. BOLD signal with N parcels as rows and
  Tmax timepoints as columns.

## Value

Numeric matrix \[N × Tmax\]. Instantaneous phases in radians.

## References

Cabral, J. et al. (2017). Cognitive performance in healthy older adults
relates to spontaneous switching between states of functional
connectivity during rest. *Scientific Reports*, 7(1), 5135.
[doi:10.1038/s41598-017-05425-7](https://doi.org/10.1038/s41598-017-05425-7)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
phases <- hilbert_phases(ts)
dim(phases)  # 10 x 200
#> [1]  10 200
```
