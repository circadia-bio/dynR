# 🧠 dynR

An R port of the Python [dynfc](https://github.com/LucasFranca/dynfc)
library for computing **dynamic connectivity (dynFC)** representations
from multivariate neurophysiological timeseries — BOLD fMRI, EEG, LFP,
and related signals.

[![R CMD
CHECK](https://github.com/circadia-bio/dynR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/circadia-bio/dynR/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://dynr.circadia-lab.uk/LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A54.1-276DC3.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.1.5-lightgrey)](https://github.com/circadia-bio/dynR)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

------------------------------------------------------------------------

> ⚠️ **dynR is in early development and has not been formally tested.**
> The API may change without notice, estimation results have not yet
> been validated against a reference implementation, and the package has
> not undergone peer review. Use with caution and verify outputs
> independently before using in any research context.

------------------------------------------------------------------------

## 📖 What is dynR?

`dynR` computes dynFC representations from preprocessed multivariate
timeseries, providing the upstream computation layer for dynamic
connectivity analysis. It is a full R port of the Python
[`dynfc`](https://github.com/LucasFranca/dynfc) library, motivated by
reproducibility: R + `renv` provides a more stable long-term environment
than Python dependency chains for research pipelines.

Although the bundled example data and several vignettes use BOLD fMRI,
all methods are applicable to any band-limited neurophysiological signal
where phase relationships or pairwise correlations carry meaningful
information — including **EEG**, **LFP**, and **MEG**.

The outputs of `dynR` feed directly into
[`stateR`](https://github.com/CoDe-Neuro/stateR) for brain state
quantification (fractional occupancy, dwell time, Markov transitions).

------------------------------------------------------------------------

## 🔁 Pipeline position

    Multivariate timeseries  [N × Tmax]
            │
            ▼
         dynR                 ← this package
      (dynFC representations)
            │
            ▼
         stateR
      (brain state metrics: FO, dwell time, Markov transitions)

------------------------------------------------------------------------

## ✨ Features

### Pipelines

| Function | Description |
|----|----|
| [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md) | Filter → Hilbert → dPL + LEiDA → Kuramoto in one call; returns `dynR_leida` |
| [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md) | Filter → sliding-window FC → cofluctuations in one call; returns `dynR_sw` |

### Phase-based methods

| Function | Description |
|----|----|
| [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md) | Instantaneous phase extraction via the analytic signal |
| [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md) | Dynamic phase-locking matrices (dPL) + LEiDA eigenvectors |
| [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md) | Leading eigenvector decomposition (LEiDA) |
| [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md) | Kuramoto order parameter, metastability, Shannon entropy |

### Correlation-based methods

| Function | Description |
|----|----|
| [`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md) | Sliding-window Pearson correlation matrices |
| [`cofluct()`](https://dynr.circadia-lab.uk/reference/cofluct.md) | Edge-centric cofluctuation time series + RSS |
| [`corr_corr()`](https://dynr.circadia-lab.uk/reference/corr_corr.md) | Correlation-of-correlations (FC recurrence) matrix |

### Utilities

| Function | Description |
|----|----|
| [`bandpass_filter()`](https://dynr.circadia-lab.uk/reference/bandpass_filter.md) | Zero-phase Butterworth bandpass filter |
| [`shannon_entropy()`](https://dynr.circadia-lab.uk/reference/shannon_entropy.md) | Shannon entropy with optional bit-depth discretisation |
| [`do_euclid()`](https://dynr.circadia-lab.uk/reference/do_euclid.md) | Euclidean distance between consecutive trajectory points |

### State dynamics

| Function | Description |
|----|----|
| [`dyn_transitions()`](https://dynr.circadia-lab.uk/reference/dyn_transitions.md) | First-order Markov transition probabilities between brain states |

### Visualisation

| Function | Description |
|----|----|
| [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md) | FC matrix heatmap (dynR diverging palette) |
| [`plot_synchrony()`](https://dynr.circadia-lab.uk/reference/plot_synchrony.md) | Kuramoto R(t) time series |
| [`plot_state_sequence()`](https://dynr.circadia-lab.uk/reference/plot_state_sequence.md) | Brain state tile plot |
| [`plot.dynR_leida()`](https://dynr.circadia-lab.uk/reference/plot.dynR_leida.md) | S3 plot method for [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md) output |
| [`plot.dynR_sw()`](https://dynr.circadia-lab.uk/reference/plot.dynR_sw.md) | S3 plot method for [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md) output |

### Multi-subject analysis

| Function | Description |
|----|----|
| [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md) | Run [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md) across a list or 3-D array of subjects |
| [`batch_sw()`](https://dynr.circadia-lab.uk/reference/batch_sw.md) | Run [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md) across a list or 3-D array of subjects |
| [`stack_leida()`](https://dynr.circadia-lab.uk/reference/stack_leida.md) | Stack LEiDA eigenvectors for cross-subject K-means clustering |
| [`stack_synchrony()`](https://dynr.circadia-lab.uk/reference/stack_synchrony.md) | Tidy long-format synchrony table across subjects |

------------------------------------------------------------------------

## ⚡ Performance

All main compute paths have compiled backends — no Python, no external
numerical libraries beyond those bundled with R:

| Function | Backend | Notes |
|----|----|----|
| [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md) | Rcpp C++ | Symmetric `cos(phi_i - phi_j)`; upper triangle only |
| [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md) | Rcpp + LAPACK `dsyev` | One shared workspace across timepoints |
| [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md) | Rcpp C++ | Direct `cos`/`sin` accumulation; no complex alloc |
| [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md) | [`mvfft()`](https://rdrr.io/r/stats/fft.html) | Two matrix FFT calls replace N per-parcel loops |
| [`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md) | Rcpp C++ | Direct Pearson; t-outer loop for column-major cache |

All backends include parity tests against their R references
(bit-perfect or \< 1e-10, depending on the algorithm).

------------------------------------------------------------------------

## 🚀 Installation

``` r

# From r-universe (recommended — no GitHub token needed)
install.packages("dynR", repos = c(
  "https://circadia-bio.r-universe.dev",
  "https://cloud.r-project.org"
))

# Or from GitHub
# install.packages("pak")
pak::pak("circadia-bio/dynR")
```

------------------------------------------------------------------------

## 📦 Quick example

``` r

library(dynR)

# Simulated timeseries: 80 channels, 300 timepoints
set.seed(42)
ts <- matrix(rnorm(80 * 300), nrow = 80)

# Pipeline wrappers — full analysis in one call
# Supply flp/fhi/delt for your modality; pass filter = FALSE if pre-filtered
res_leida <- leida_pipeline(ts, flp = 0.01, fhi = 0.1, delt = 2)  # fMRI example
res_sw    <- sw_pipeline(ts, window = 30, step = 5, flp = 0.01, fhi = 0.1, delt = 2)

# Inspect
res_leida            # <dynR_leida> — N, Tmax, metastability, entropy
plot(res_leida)      # Kuramoto R(t)
plot(res_sw, "fc")   # Mean sliding-window FC

# Cluster LEiDA eigenvectors and visualise the state sequence
km     <- kmeans(res_leida$leida, centers = 5, nstart = 100)
plot_state_sequence(km$cluster)

# Multi-subject: list or [N × Tmax × subjects] array
batch    <- batch_leida(list(sub01 = ts, sub02 = ts), filter = FALSE)
df_sync  <- stack_synchrony(batch)   # tidy: subject / timepoint / synchrony
leida_all <- stack_leida(batch)      # ready for cross-subject kmeans
```

Step-by-step (without pipeline wrappers)

``` r

# 1. Bandpass filter (fMRI: TR = 2 s, 0.01–0.1 Hz)
ts_filt <- t(apply(ts, 1, bandpass_filter, flp = 0.01, fhi = 0.1, delt = 2))

# 2. Phase-based: LEiDA + Kuramoto
phases <- hilbert_phases(ts_filt)
dpl    <- dyn_phase_lock(phases)   # dpl$leida: [Tmax-20 × 80]
kop    <- kuramoto(phases)         # kop$metastability, kop$entropy

# 3. Correlation-based: sliding-window FC + cofluctuations
sw <- corr_slide(ts_filt, window = 30, step = 5)
ec <- cofluct(ts_filt)             # ec$edge_ts, ec$rss
```

------------------------------------------------------------------------

## 📐 Data conventions

All timeseries inputs follow:

> **rows = channels/parcels (N), columns = timepoints (Tmax)**

This matches the `[N, Tmax]` convention of the original Python `dynfc`
package.

------------------------------------------------------------------------

## 👥 Authors

| Role | Name | Affiliation |
|----|----|----|
| Author, maintainer | Lucas G. S. França | Northumbria University / Circadia Lab |
| Author | Mario Leocadio-Miguel | Northumbria University / Circadia Lab |
| Author | Dafnis Batallé | King’s College London / CoDe-Neuro Lab |

------------------------------------------------------------------------

## 🤝 Related tools

**Circadia Lab ecosystem:** - 🕐
[**zeitR**](https://github.com/circadia-bio/zeitR) — wrist actigraphy
analysis - 😴 [**hypnoR**](https://github.com/circadia-bio/hypnoR) —
hypnogram handling and sleep architecture - 🔗
[**syncR**](https://github.com/circadia-bio/syncR) — ecosystem
integrator

**CoDe-Neuro ecosystem:** - 🧪
[**stateR**](https://github.com/CoDe-Neuro/stateR) — brain state metrics
(FO, dwell time, Markov) — consumes `dynR` output - 🔬
[**ptestR**](https://github.com/CoDe-Neuro/ptestR) — permutation-based
significance testing

------------------------------------------------------------------------

## 📄 Licence

MIT © 2026 Lucas G. S. França, Mario Leocadio-Miguel, Dafnis Batallé
