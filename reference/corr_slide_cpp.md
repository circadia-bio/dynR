# Sliding-window Pearson correlation (C++ backend)

Internal workhorse called by
[`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md).
Computes Pearson correlation matrices directly on the input matrix
without allocating per-window submatrices or calling R's
[`cor()`](https://rdrr.io/r/stats/cor.html).

## Usage

``` r
corr_slide_cpp(timeseries, window, step)
```

## Arguments

- timeseries:

  NumericMatrix \[N x Tmax\].

- window:

  Integer window size in timepoints.

- step:

  Integer step between window onsets.

## Value

List with:

- corr_mats:

  NumericVector with dim \[N, N, n_windows\].

- idx:

  IntegerVector of 1-indexed window onset positions.

## Details

Loop ordering is chosen for cache efficiency on the column-major R
matrix: the innermost loop over `t` (timepoints) within a window
iterates over consecutive memory locations in each column of
`timeseries`, so the cross-product accumulation step is sequential in
memory.
