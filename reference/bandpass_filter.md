# Butterworth bandpass filter

Design and apply a zero-phase Butterworth bandpass filter to a signal.
Wraps [`gsignal::butter()`](https://rdrr.io/pkg/gsignal/man/butter.html)
and
[`gsignal::filtfilt()`](https://rdrr.io/pkg/gsignal/man/filtfilt.html)
for convenience.

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

## Examples

``` r
set.seed(1)
x <- rnorm(200)
x_filt <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 2, order = 2)
```
