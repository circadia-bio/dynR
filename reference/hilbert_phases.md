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

  Numeric matrix \[N x Tmax\]. BOLD signal with N parcels as rows and
  Tmax timepoints as columns.

## Value

Numeric matrix \[N x Tmax\]. Instantaneous phases in radians.

## Details

All N parcels are processed in two vectorised calls to
[`mvfft()`](https://rdrr.io/r/stats/fft.html) rather than an N-iteration
R loop over [`fft()`](https://rdrr.io/r/stats/fft.html): the timeseries
matrix is transposed to a Tmax x N layout so `mvfft` applies the FFT to
each parcel column simultaneously, the Hilbert multiplier is broadcast
across columns, and the inverse FFT recovers the analytic signal for all
parcels at once.

`.hilbert_r()` is retained as an internal single-vector reference
(scipy-parity validated).

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
