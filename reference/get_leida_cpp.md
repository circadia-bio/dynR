# Leading eigenvector decomposition (LEiDA) – C++ / LAPACK backend

Internal workhorse called by
[`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md).
For each phase-locking matrix slice in `sync_conn`, calls LAPACK `dsyev`
(symmetric eigendecomposition) and returns only the leading eigenvector.

## Usage

``` r
get_leida_cpp(sync_conn)
```

## Arguments

- sync_conn:

  NumericVector with `dim` attribute \[N, N, t_points\].

## Value

NumericMatrix \[t_points x N\]. Leading eigenvectors, one per row.

## Details

Differences from the R implementation:

- `dsyev` returns eigenvalues in **ascending** order, so the leading
  eigenvector is the **last** column of the output (opposite of R's
  [`eigen()`](https://rdrr.io/r/base/eigen.html) which is descending).

- The LAPACK workspace is allocated once and reused across all
  timepoints.

- Sign convention preserved: row sum forced to be non-positive.
