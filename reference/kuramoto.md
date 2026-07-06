# Kuramoto order parameter and metastability

Compute the global Kuramoto order parameter time series, the
metastability index (standard deviation of the order parameter), and
Shannon entropy of synchrony from parcel-level instantaneous phase time
series. The first and last 10 timepoints are discarded (matching
[`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md)).

## Usage

``` r
kuramoto(phases, base = 2, n_bits = 8L)
```

## Arguments

- phases:

  Numeric matrix \[N x Tmax\]. Instantaneous phases in radians, as
  returned by
  [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md).

- base:

  Numeric. Logarithm base for Shannon entropy. Default is 2.

- n_bits:

  Integer or `NULL`. Bit depth for discretising the synchrony series
  before entropy estimation. Default is 8.

## Value

A list with:

- metastability:

  Numeric scalar. Standard deviation of the Kuramoto order parameter.

- synchrony:

  Numeric vector \[Tmax-20\]. Kuramoto order parameter at each (trimmed)
  timepoint.

- entropy:

  Numeric scalar. Shannon entropy of the synchrony series.

## Details

The synchrony time series is computed by a compiled C++ backend
([`kuramoto_sync_cpp()`](https://dynr.circadia-lab.uk/reference/kuramoto_sync_cpp.md))
that accumulates cos and sin components directly, avoiding R
complex-number allocation and `vapply` dispatch overhead.

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
phases <- hilbert_phases(ts)
res <- kuramoto(phases)
res$metastability
#> [1] 0.1427287
```
