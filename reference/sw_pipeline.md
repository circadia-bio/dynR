# Sliding-window pipeline

Convenience wrapper running the full correlation-based dynamic FC
pipeline: optional Butterworth bandpass filter -\> sliding-window
Pearson correlation matrices -\> edge-centric cofluctuations. Returns a
structured `dynR_sw` object with
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Usage

``` r
sw_pipeline(
  timeseries,
  window,
  step = NULL,
  flp = NULL,
  fhi = NULL,
  delt = NULL,
  order = 2L,
  filter = TRUE
)
```

## Arguments

- timeseries:

  Numeric matrix \[N x Tmax\].

- window:

  Integer. Window size in timepoints.

- step:

  Integer. Step between window onsets. Default: `window`
  (non-overlapping windows).

- flp:

  Numeric. Low-pass cutoff (Hz). No default — required when
  `filter = TRUE`.

- fhi:

  Numeric. High-pass cutoff (Hz). No default — required when
  `filter = TRUE`.

- delt:

  Numeric. Sampling interval in seconds. No default — required when
  `filter = TRUE`.

- order:

  Integer. Butterworth filter order. Default `2L`.

- filter:

  Logical. Apply bandpass filter? Default `TRUE`.

## Value

An object of class `dynR_sw` (a named list) with elements:

- corr_mats:

  Array \[N, N, n_windows\]. Sliding FC matrices.

- idx:

  Integer vector. 1-indexed window onset positions.

- edge_ts:

  Matrix \[n_edges, Tmax\]. Edge time series.

- rss:

  Numeric vector \[Tmax\]. Root-sum-square cofluctuation.

- window:

  Integer. Window size used.

- step:

  Integer. Step size used.

- N:

  Integer. Number of channels/parcels.

- Tmax:

  Integer. Number of timepoints.

## Details

dynR is modality-agnostic. `flp`, `fhi`, and `delt` have no defaults and
**must** be supplied when `filter = TRUE`. Pass `filter = FALSE` if the
timeseries is already band-limited.

## See also

[`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md),
[`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md),
[`cofluct()`](https://dynr.circadia-lab.uk/reference/cofluct.md),
[`plot.dynR_sw()`](https://dynr.circadia-lab.uk/reference/plot.dynR_sw.md)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10)

# filter = FALSE: timeseries already band-limited
res <- sw_pipeline(ts, window = 20, filter = FALSE)
res
#> <dynR_sw>
#>   Parcels:     10 
#>   Timepoints:  200 
#>   Window:      20 timepoints
#>   Step:        20 timepoints
#>   N windows:   10 
#>   N edges:     45 

# fMRI BOLD (TR = 2 s)
if (FALSE) { # \dontrun{
res <- sw_pipeline(ts, window = 20, flp = 0.01, fhi = 0.1, delt = 2)
} # }
```
