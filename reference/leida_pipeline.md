# LEiDA pipeline

Convenience wrapper running the full phase-based dynamic FC pipeline:
optional Butterworth bandpass filter -\> Hilbert transform -\> dynamic
phase-locking matrices + LEiDA eigenvectors -\> Kuramoto order
parameter. Returns a structured `dynR_leida` object with
[`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Usage

``` r
leida_pipeline(
  timeseries,
  flp = NULL,
  fhi = NULL,
  delt = NULL,
  order = 2L,
  filter = TRUE
)
```

## Arguments

- timeseries:

  Numeric matrix \[N x Tmax\]. Rows = channels/parcels, columns =
  timepoints.

- flp:

  Numeric. Low-pass cutoff frequency (Hz). No default — must be set
  explicitly when `filter = TRUE`.

- fhi:

  Numeric. High-pass cutoff frequency (Hz). No default — must be set
  explicitly when `filter = TRUE`.

- delt:

  Numeric. Sampling interval in seconds (`1 / sampling_rate`). No
  default — must be set explicitly when `filter = TRUE`.

- order:

  Integer. Butterworth filter order. Default `2L`.

- filter:

  Logical. Apply bandpass filter before phase extraction? Default
  `TRUE`. Set to `FALSE` if the timeseries is already filtered.

## Value

An object of class `dynR_leida` (a named list) with elements:

- leida:

  Matrix \[Tmax-20, N\]. LEiDA eigenvectors.

- sync_conn:

  Array \[N, N, Tmax-20\]. Instantaneous phase-locking matrices.

- phases:

  Matrix \[N, Tmax\]. Instantaneous phases in radians.

- synchrony:

  Numeric vector \[Tmax-20\]. Kuramoto R(t).

- metastability:

  Numeric. Standard deviation of synchrony.

- entropy:

  Numeric. Shannon entropy of synchrony series (bits).

- N:

  Integer. Number of channels/parcels.

- Tmax:

  Integer. Number of timepoints.

## Details

dynR is modality-agnostic. `flp`, `fhi`, and `delt` have no defaults and
**must** be supplied when `filter = TRUE`. Set them for your recording
modality (see examples). Pass `filter = FALSE` if the timeseries is
already band-limited.

## See also

[`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md),
[`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md),
[`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md),
[`plot.dynR_leida()`](https://dynr.circadia-lab.uk/reference/plot.dynR_leida.md)

## Examples

``` r
set.seed(1)
ts <- matrix(rnorm(10 * 200), nrow = 10)

# filter = FALSE: timeseries already band-limited
res <- leida_pipeline(ts, filter = FALSE)
res
#> <dynR_leida>
#>   Parcels:       10 
#>   Timepoints:    200 
#>   LEiDA frames:  180 
#>   Metastability: 0.1427 
#>   Entropy:       6.4623 bits

# fMRI BOLD (TR = 2 s)
if (FALSE) { # \dontrun{
res <- leida_pipeline(ts, flp = 0.01, fhi = 0.1, delt = 2)
} # }

# EEG alpha band (250 Hz)
if (FALSE) { # \dontrun{
res <- leida_pipeline(ts, flp = 8, fhi = 12, delt = 1 / 250)
} # }
```
