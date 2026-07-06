# Kuramoto order parameter time series (C++ backend)

Internal workhorse called by
[`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md). For
each trimmed timepoint, computes the global Kuramoto order parameter
R(t):

## Usage

``` r
kuramoto_sync_cpp(phases)
```

## Arguments

- phases:

  NumericMatrix \[N x Tmax\]. Instantaneous phases in radians.

## Value

NumericVector \[Tmax-20\]. Kuramoto order parameter per timepoint.

## Details

\$\$R(t) = \frac{1}{N} \left\| \sum\_{j=1}^{N} e^{i\phi_j(t)}
\right\|\$\$

The complex exponential is decomposed as \\\cos\phi + i\sin\phi\\; the
modulus is \\\sqrt{(\sum\cos)^2 + (\sum\sin)^2} / N\\, avoiding any R
complex-number allocation.

The first and last 10 timepoints are discarded (matching
[`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md)).
