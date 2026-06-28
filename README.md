# 🧠 dynR

An R package for **dynamic functional connectivity (dynFC)** analysis of BOLD fMRI timeseries — sliding-window correlations, edge-centric cofluctuations, instantaneous phase-locking, LEiDA leading eigenvectors, and the Kuramoto order parameter.

[![R CMD CHECK](https://github.com/CoDe-Neuro/dynR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CoDe-Neuro/dynR/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)
[![R](https://img.shields.io/badge/R-≥4.1-276DC3.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/version-0.1.0-lightgrey)](https://github.com/CoDe-Neuro/dynR)

---

## 📖 What is dynR?

`dynR` is the upstream dynamic FC computation layer in the CoDe-Neuro R ecosystem.
It picks up preprocessed BOLD timeseries — typically from
[fMRIPrep](https://fmriprep.org) derivatives, or the output of
[`boldR`](https://github.com/circadia-bio/boldR) — and transforms them into
dynFC representations ready for brain state analysis in
[`stateR`](https://github.com/CoDe-Neuro/stateR).

`dynR` is a full R port of the Python [`dynfc`](https://github.com/LucasFranca/dynfc)
package, motivated by reproducibility: R + `renv` provides a more stable
long-term environment than Python dependency chains for research pipelines.

---

## 🔁 Pipeline position

```
BOLD timeseries
      │
      ▼
   dynR                 ← this package
  (dynFC representations)
      │
      ▼
   stateR
  (brain state metrics: FO, dwell time, Markov transitions)
```

`dynR` outputs are designed to feed directly into `stateR::nest_fo()`,
`stateR::nest_dwell()`, and `stateR::clusters_markov()`.

---

## ✨ Features

### Phase-based methods
| Function | Description |
|---|---|
| `hilbert_phases()` | Instantaneous phase extraction via the analytic signal |
| `dyn_phase_lock()` | Dynamic phase-locking matrices (dPL) + LEiDA eigenvectors |
| `get_leida()` | Leading eigenvector decomposition (LEiDA) |
| `kuramoto()` | Kuramoto order parameter, metastability, Shannon entropy |

### Correlation-based methods
| Function | Description |
|---|---|
| `corr_slide()` | Sliding-window Pearson correlation matrices |
| `cofluct()` | Edge-centric cofluctuation time series + RSS |
| `corr_corr()` | Correlation-of-correlations matrix |

### Utilities
| Function | Description |
|---|---|
| `bandpass_filter()` | Zero-phase Butterworth bandpass filter |
| `shannon_entropy()` | Shannon entropy with optional bit-depth discretisation |
| `do_euclid()` | Euclidean distance between consecutive trajectory points |

---

## 🚀 Installation

```r
# install.packages("pak")
pak::pak("CoDe-Neuro/dynR")
```

---

## 📦 Quick example

```r
library(dynR)

# Simulated BOLD: 80 parcels, 300 timepoints
set.seed(42)
ts <- matrix(rnorm(80 * 300), nrow = 80, ncol = 300)

# 1. Bandpass filter (TR = 2 s)
ts_filt <- apply(ts, 1, bandpass_filter, flp = 0.01, fhi = 0.1, delt = 2)
ts_filt <- t(ts_filt)

# 2. Extract instantaneous phases
phases <- hilbert_phases(ts_filt)

# 3. Dynamic phase-locking + LEiDA eigenvectors
dpl <- dyn_phase_lock(phases)
# dpl$sync_conn  — [80, 80, 280] array of phase-locking matrices
# dpl$leida      — [280, 80] matrix of leading eigenvectors

# 4. Cluster LEiDA vectors (→ feed into stateR)
km <- kmeans(dpl$leida, centers = 5, nstart = 100)
# km$cluster is the state sequence for stateR

# 5. Kuramoto order parameter
kop <- kuramoto(phases)
kop$metastability
kop$entropy

# 6. Sliding-window FC
sw <- corr_slide(ts_filt, window = 30, step = 5)
dim(sw$corr_mats)  # 80 x 80 x n_windows

# 7. Edge-centric cofluctuations
ec <- cofluct(ts_filt)
# ec$edge_ts — edge time series [n_edges x 300]
# ec$rss     — root-sum-square cofluctuation [300]
```

---

## 🔧 Function reference

### `hilbert_phases(timeseries)`

Computes instantaneous phases via the Hilbert transform.
Input: `[N × Tmax]` matrix. Output: `[N × Tmax]` phase matrix in radians.
Each parcel is demeaned before transformation.

### `dyn_phase_lock(phases)`

Computes instantaneous phase-locking matrices from phase time series.
Trims 10 timepoints from each end (edge effects).
Returns `sync_conn` `[N, N, Tmax-20]` and `leida` `[Tmax-20, N]`.

### `get_leida(sync_conn)`

Extracts the leading eigenvector from each phase-locking matrix.
Sign-normalised so that the sum of each eigenvector is ≤ 0.
Output: `[Tmax × N]`.

### `kuramoto(phases, base = 2, n_bits = 8)`

Global Kuramoto order parameter time series (trimmed), metastability index
(standard deviation), and Shannon entropy of synchrony.

### `corr_slide(timeseries, window, step = NULL)`

Pearson correlation matrices over sliding windows. Non-overlapping by default
(`step = window`). Returns `corr_mats` `[N, N, n_windows]` and onset `idx`.

### `cofluct(timeseries, k = 1)`

Edge time series (element-wise product of z-scored parcel pairs) and RSS
cofluctuation vector. Returns `edge_ts` `[n_edges × Tmax]` and `rss` `[Tmax]`.

### `corr_corr(timeseries, k = 1)`

Correlation of correlations: `[Tmax × Tmax]` matrix of pairwise correlations
between edge time series.

### `bandpass_filter(x, flp, fhi, delt, order = 2)`

Zero-phase Butterworth bandpass filter via `gsignal::butter()` +
`gsignal::filtfilt()`.

### `shannon_entropy(x, base = 2, n_bits = NULL)`

Shannon entropy from a numeric vector. Optional bit-depth discretisation.

### `do_euclid(x)`

Row-wise Euclidean distance between consecutive points in a matrix.

---

## 📐 Data conventions

All timeseries inputs follow the convention:

> **rows = parcels/voxels (N), columns = timepoints (Tmax)**

This matches the `[N, Tmax]` convention used in the original Python `dynfc` package.

---

## 👥 Authors

| Role | Name | Affiliation |
|---|---|---|
| Author, maintainer | Lucas G. S. França | Northumbria University / CoDe-Neuro Lab |
| Author | Dafnis Batallé | King's College London / CoDe-Neuro Lab |

Part of the [CoDe-Neuro Lab](https://github.com/CoDe-Neuro) at King's College London.

---

## 🤝 Related tools

- 🧪 [**stateR**](https://github.com/CoDe-Neuro/stateR) — Brain state metrics (FO, dwell time, Markov) — consumes `dynR` output
- 🔬 [**ptestR**](https://github.com/CoDe-Neuro/ptestR) — Permutation-based significance testing for `stateR` outputs
- 🧲 [**boldR**](https://github.com/circadia-bio/boldR) — fMRI BOLD preprocessing and parcellation; feeds `dynR`

---

## 📄 Licence

MIT © 2026 Lucas G. S. França, Dafnis Batallé
