# Getting started with dynR

## What is dynamic functional connectivity?

Resting-state fMRI measures the BOLD signal across brain regions over
time. **Static** functional connectivity (FC) summarises this as a
single correlation matrix — the average co-activation pattern across the
entire scan. **Dynamic** functional connectivity (dynFC) asks a
different question: does that pattern change over time, and if so, how?

`dynR` computes dynFC representations from preprocessed BOLD timeseries.
It implements two complementary families of methods:

| Family | Key functions | What it captures |
|----|----|----|
| Phase-based | [`hilbert_phases()`](https://CoDe-Neuro.github.io/dynR/reference/hilbert_phases.md), [`dyn_phase_lock()`](https://CoDe-Neuro.github.io/dynR/reference/dyn_phase_lock.md), [`kuramoto()`](https://CoDe-Neuro.github.io/dynR/reference/kuramoto.md) | Instantaneous synchrony and phase relationships |
| Correlation-based | [`corr_slide()`](https://CoDe-Neuro.github.io/dynR/reference/corr_slide.md), [`cofluct()`](https://CoDe-Neuro.github.io/dynR/reference/cofluct.md), [`corr_corr()`](https://CoDe-Neuro.github.io/dynR/reference/corr_corr.md) | Windowed or edge-level co-activation |

The outputs of `dynR` feed directly into
[`stateR`](https://github.com/CoDe-Neuro/stateR) for brain state
quantification (fractional occupancy, dwell time, Markov transitions).

------------------------------------------------------------------------

## Package data

`dynR` ships with a resting-state BOLD timeseries from the
[edge-ts](https://github.com/brain-networks/edge-ts) repository,
parcellated into 200 regions:

``` r

data(ts)   # [200 x 600] BOLD timeseries
data(fc)   # [200 x 200] static FC matrix (ground truth)

dim(ts)
#> [1] 200 600
dim(fc)
#> [1] 200 200
```

200 parcels, 600 timepoints. We will use this throughout the vignettes.

------------------------------------------------------------------------

## Pipeline overview

A typical dynR workflow looks like this:

    BOLD timeseries [N × Tmax]
            │
            ├─── bandpass_filter()       ← optional preprocessing
            │
            ├─── Phase-based ────────────────────────────────────────
            │     hilbert_phases()       → instantaneous phases [N × Tmax]
            │     dyn_phase_lock()       → phase-locking matrices + LEiDA vectors
            │     kuramoto()             → synchrony, metastability, entropy
            │
            └─── Correlation-based ──────────────────────────────────
                  corr_slide()          → sliding-window FC [N × N × windows]
                  cofluct()             → edge time series + RSS
                  corr_corr()           → correlation-of-correlations [Tmax × Tmax]
                        │
                        ▼
                  stateR (K-means → brain states)

------------------------------------------------------------------------

## Quick example: phase-based pipeline

``` r

# 1. Extract instantaneous phases via the Hilbert transform
phases <- hilbert_phases(ts)
dim(phases)  # same as ts: [200 x 600]
#> [1] 200 600

# 2. Dynamic phase-locking matrices + leading eigenvectors (LEiDA)
dpl <- dyn_phase_lock(phases)
dim(dpl$sync_conn)  # [200 x 200 x 580]  (trimmed 10 timepoints each end)
#> [1] 200 200 580
dim(dpl$leida)      # [580 x 200]         leading eigenvectors
#> [1] 580 200

# 3. Kuramoto order parameter
kop <- kuramoto(phases)
cat("Metastability:", round(kop$metastability, 4), "\n")
#> Metastability: 0.0796
cat("Entropy:      ", round(kop$entropy, 4), "\n")
#> Entropy:       6.1484
```

The `leida` matrix — one leading eigenvector per timepoint — is ready to
pass to K-means for brain state discovery.

------------------------------------------------------------------------

## Quick example: sliding-window pipeline

``` r

# Sliding-window correlation: 30-timepoint windows, step of 5
sw <- corr_slide(ts, window = 30, step = 5)
cat("Number of windows:", dim(sw$corr_mats)[3], "\n")
#> Number of windows: 115

# Edge-centric cofluctuations
ec <- cofluct(ts)
cat("Edge time series shape:", nrow(ec$edge_ts), "x", ncol(ec$edge_ts), "\n")
#> Edge time series shape: 19900 x 600
cat("RSS range: [", round(min(ec$rss), 2), ",", round(max(ec$rss), 2), "]\n")
#> RSS range: [ 12.4 , 499.65 ]
```

------------------------------------------------------------------------

## Validating against static FC

A useful sanity check: a single window spanning the full timeseries
should reproduce the static FC matrix exactly.

``` r

full_window <- corr_slide(ts, window = ncol(ts))
max_diff    <- max(abs(full_window$corr_mats[, , 1] - fc))
cat("Max deviation from static FC:", formatC(max_diff, format = "e"), "\n")
#> Max deviation from static FC: 8.8818e-16
```

------------------------------------------------------------------------

## What comes next

- [`vignette("phase-based-fc")`](https://CoDe-Neuro.github.io/dynR/articles/phase-based-fc.md)
  — Hilbert transform, LEiDA, and the Kuramoto order parameter in depth
- [`vignette("sliding-window-fc")`](https://CoDe-Neuro.github.io/dynR/articles/sliding-window-fc.md)
  — sliding-window correlation and edge-centric cofluctuations
- Passing `dpl$leida` to
  [`kmeans()`](https://rdrr.io/r/stats/kmeans.html) and `stateR` for
  brain state quantification
