# Instantaneous phase-locking matrices (C++ backend)

Internal workhorse called by
[`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md).
Iterates over the trimmed timepoints (dropping the first and last 10)
and fills an \[N, N, n_T\] array with \\\cos(\phi_i - \phi_j)\\ for
every parcel pair.

## Usage

``` r
dyn_phase_lock_cpp(phases)
```

## Arguments

- phases:

  NumericMatrix \[N x Tmax\]. Instantaneous phases in radians, as
  returned by
  [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md).

## Value

A NumericVector with a `dim` attribute \[N, N, n_T\] (i.e. a 3-D array),
where `n_T = Tmax - 20`.

## Details

Exploits the symmetry \\\cos(a-b) = \cos(b-a)\\: only the upper triangle
is computed via `std::cos`; the diagonal is set to 1 directly, and
values are mirrored to the lower triangle. This halves the number of
trigonometric evaluations relative to a naive double loop.
