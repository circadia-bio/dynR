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

Two internal code paths are used depending on `order`:

- **`order <= 3`** (typical fMRI use): direct-form b/a coefficients with
  the companion-matrix initial-condition method. scipy `filtfilt`-parity
  validated to `< 1e-9` for standard fMRI parameters.

- **`order >= 4`**: second-order sections (SOS). The companion matrix
  for the full-order b/a representation becomes ill-conditioned for
  high-order filters with poles close to z = 1 (very low lower cutoff).
  SOS avoids this by decomposing the filter into biquad sections, each
  with a 2x2 companion matrix that is always well-conditioned.

Filter coefficients are computed by `.butter_bandpass_zpk()`, a bespoke
port of `scipy.signal.butter()` that requires no external packages.

## Examples

``` r
set.seed(1)
x <- rnorm(200)
x_filt <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 2, order = 2)
```
