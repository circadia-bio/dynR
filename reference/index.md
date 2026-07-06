# Package index

## Data

Example datasets bundled with dynR.

- [`ts`](https://dynr.circadia-lab.uk/reference/ts.md) : BOLD fMRI
  timeseries (200 parcels, 600 timepoints)
- [`fc`](https://dynr.circadia-lab.uk/reference/fc.md) : Functional
  connectivity matrix (200 parcels)

## Pipelines

High-level wrappers running the full phase-based or correlation-based
pipeline in one call.

- [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md)
  : LEiDA pipeline
- [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md)
  : Sliding-window pipeline

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

## Visualisation

Plot FC matrices, synchrony time series, and brain state sequences.

- [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md) :
  Plot a functional connectivity matrix
- [`plot_synchrony()`](https://dynr.circadia-lab.uk/reference/plot_synchrony.md)
  : Plot Kuramoto synchrony time series
- [`plot_state_sequence()`](https://dynr.circadia-lab.uk/reference/plot_state_sequence.md)
  : Plot a brain state sequence
- [`plot(`*`<dynR_leida>`*`)`](https://dynr.circadia-lab.uk/reference/plot.dynR_leida.md)
  : Plot method for dynR_leida objects
- [`plot(`*`<dynR_sw>`*`)`](https://dynr.circadia-lab.uk/reference/plot.dynR_sw.md)
  : Plot method for dynR_sw objects

## Multi-subject analysis

Batch processing across participants and stacking for cross-subject
clustering.

- [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md)
  : Batch LEiDA pipeline across subjects
- [`batch_sw()`](https://dynr.circadia-lab.uk/reference/batch_sw.md) :
  Batch sliding-window pipeline across subjects
- [`stack_leida()`](https://dynr.circadia-lab.uk/reference/stack_leida.md)
  : Stack LEiDA eigenvectors across subjects
- [`stack_synchrony()`](https://dynr.circadia-lab.uk/reference/stack_synchrony.md)
  : Stack synchrony vectors across subjects

## Utilities

Entropy and distance helpers.

- [`shannon_entropy()`](https://dynr.circadia-lab.uk/reference/shannon_entropy.md)
  : Shannon entropy
- [`do_euclid()`](https://dynr.circadia-lab.uk/reference/do_euclid.md) :
  Euclidean distance between consecutive points
