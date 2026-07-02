# Butterworth bandpass filter

Design and apply a zero-phase Butterworth bandpass filter to a signal,
using scipy-compatible steady-state initial conditions (equivalent to
`scipy.signal.sosfiltfilt` / `scipy.signal.filtfilt` with
`padtype = "odd"`).

## Usage

``` r
bandpass_filter(x, flp, fhi, delt, order = 2)
```

## Arguments

- x:

  Numeric vector. Signal to be filtered.

- flp:

  Numeric. Low-pass cutoff frequency (Hz).

- fhi:

  Numeric. High-pass cutoff frequency (Hz).

- delt:

  Numeric. Sampling interval in seconds (i.e. TR for fMRI).

- order:

  Integer. Filter order. Default is 2.

## Value

Numeric vector. Zero-phase filtered signal, same length as `x`.

## Details

[`gsignal::filtfilt()`](https://rdrr.io/pkg/gsignal/man/filtfilt.html)
uses zero initial conditions, which produces edge transients up to ~0.24
signal units on real fMRI data. This implementation replicates scipy's
`lfilter_zi` approach: both the forward and backward passes are
initialised at steady state scaled by the first sample of each pass,
reducing edge error to machine precision.

## Examples

``` r
set.seed(1)
x <- rnorm(200)
x_filt <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 2, order = 2)
```
