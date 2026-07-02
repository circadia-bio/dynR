# Package index

## Data

Example datasets bundled with dynR.

- [`ts`](https://dynr.circadia-lab.uk/reference/ts.md) : BOLD fMRI
  timeseries (200 parcels, 600 timepoints)
- [`fc`](https://dynr.circadia-lab.uk/reference/fc.md) : Functional
  connectivity matrix (200 parcels)

## Signal preprocessing

Filter and prepare BOLD timeseries prior to dynFC analysis.

- [`bandpass_filter()`](https://dynr.circadia-lab.uk/reference/bandpass_filter.md)
  : Butterworth bandpass filter

## Phase-based methods

Instantaneous phase extraction, phase-locking, and Kuramoto order
parameter.

- [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md)
  : Hilbert transform phase extraction
- [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md)
  : Dynamic phase-locking matrix (dPL)
- [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md) :
  Leading eigenvector decomposition (LEiDA)
- [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md) :
  Kuramoto order parameter and metastability

## Correlation-based methods

Sliding-window FC and edge-centric cofluctuation.

- [`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md)
  : Sliding window correlation
- [`cofluct()`](https://dynr.circadia-lab.uk/reference/cofluct.md) :
  Edge-centric cofluctuation analysis
- [`corr_corr()`](https://dynr.circadia-lab.uk/reference/corr_corr.md) :
  Correlation of correlations matrix

## State dynamics

Quantify the temporal structure of brain-state sequences.

- [`dyn_transitions()`](https://dynr.circadia-lab.uk/reference/dyn_transitions.md)
  : State transition probabilities (Markov analysis)

## Utilities

Entropy and distance helpers.

- [`shannon_entropy()`](https://dynr.circadia-lab.uk/reference/shannon_entropy.md)
  : Shannon entropy
- [`do_euclid()`](https://dynr.circadia-lab.uk/reference/do_euclid.md) :
  Euclidean distance between consecutive points
