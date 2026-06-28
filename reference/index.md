# Package index

## Data

Example datasets bundled with dynR.

- [`ts`](https://CoDe-Neuro.github.io/dynR/reference/ts.md) : BOLD fMRI
  timeseries (200 parcels, 600 timepoints)
- [`fc`](https://CoDe-Neuro.github.io/dynR/reference/fc.md) : Functional
  connectivity matrix (200 parcels)

## Signal preprocessing

Filter and prepare BOLD timeseries prior to dynFC analysis.

- [`bandpass_filter()`](https://CoDe-Neuro.github.io/dynR/reference/bandpass_filter.md)
  : Butterworth bandpass filter

## Phase-based methods

Instantaneous phase extraction, phase-locking, and Kuramoto order
parameter.

- [`hilbert_phases()`](https://CoDe-Neuro.github.io/dynR/reference/hilbert_phases.md)
  : Hilbert transform phase extraction
- [`dyn_phase_lock()`](https://CoDe-Neuro.github.io/dynR/reference/dyn_phase_lock.md)
  : Dynamic phase-locking matrix (dPL)
- [`get_leida()`](https://CoDe-Neuro.github.io/dynR/reference/get_leida.md)
  : Leading eigenvector decomposition (LEiDA)
- [`kuramoto()`](https://CoDe-Neuro.github.io/dynR/reference/kuramoto.md)
  : Kuramoto order parameter and metastability

## Correlation-based methods

Sliding-window FC and edge-centric cofluctuation.

- [`corr_slide()`](https://CoDe-Neuro.github.io/dynR/reference/corr_slide.md)
  : Sliding window correlation
- [`cofluct()`](https://CoDe-Neuro.github.io/dynR/reference/cofluct.md)
  : Edge-centric cofluctuation analysis
- [`corr_corr()`](https://CoDe-Neuro.github.io/dynR/reference/corr_corr.md)
  : Correlation of correlations matrix

## Utilities

Entropy and distance helpers.

- [`shannon_entropy()`](https://CoDe-Neuro.github.io/dynR/reference/shannon_entropy.md)
  : Shannon entropy
- [`do_euclid()`](https://CoDe-Neuro.github.io/dynR/reference/do_euclid.md)
  : Euclidean distance between consecutive points
