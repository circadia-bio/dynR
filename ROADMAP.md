# dynR roadmap

dynR is modality-agnostic — BOLD fMRI, EEG, LFP, MEG — and the roadmap
reflects that. Upcoming milestones add the connectivity measures and
analysis patterns that non-fMRI recordings actually require, rather than
simply extending fMRI-centric defaults to other modalities.

------------------------------------------------------------------------

## v0.2.0 — cross-modal connectivity

Adds EEG/LFP/MEG-appropriate connectivity measures that operate directly
on the phase and amplitude outputs already produced by
[`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md).
No new dependencies; all three get Rcpp backends.

| Function | Description |
|----|----|
| `compute_pli()` | Phase lag index — discards zero-lag synchrony via `sign(sin(Δφ))` |
| `compute_wpli()` | Weighted PLI — standard EEG connectivity; de Haan et al. 2011 |
| `amplitude_envelope()` | Expose Hilbert amplitudes (already computed, currently discarded) |
| `aec()` | Amplitude envelope correlation — MEG / source-space standard |

**Why PLI/wPLI over PLV for EEG?**
[`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md)
computes the phase locking value (`cos(Δφ)`), which is contaminated by
volume conduction in scalp EEG — two electrodes picking up the same
source appear perfectly synchronised at zero lag. PLI and wPLI discard
or downweight zero-lag contributions, making them scientifically
defensible for EEG without requiring source localisation.

**Why AEC?** Phase-based measures assume the signal is oscillatory. For
broadband LFP or source-reconstructed MEG, the correlation of Hilbert
amplitude envelopes (AEC) is often preferred and is what
`amplitude_envelope()` enables.

------------------------------------------------------------------------

## v0.2.x — multi-band analysis

fMRI is effectively single-band (0.01–0.1 Hz). EEG and LFP are
multi-band — delta, theta, alpha, beta, gamma each carry distinct
information. This milestone adds wrappers that loop over named band
definitions and return structured per-band outputs.

| Function | Description |
|----|----|
| `multi_band_leida()` | [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md) across a named list of `list(flp=, fhi=)` band definitions |
| `multi_band_sw()` | [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md) across the same band list |
| Band presets | `dynr_bands_eeg()` returning delta/theta/alpha/beta/gamma defaults |

Per-band `dynR_leida` / `dynR_sw` outputs are designed to be stackable
via the existing
[`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md)
/
[`stack_leida()`](https://dynr.circadia-lab.uk/reference/stack_leida.md)
infrastructure.

------------------------------------------------------------------------

## v0.3.0 — extended methods

| Item | Notes |
|----|----|
| Sleep cycle identification | Shared upstream concern with mrpheus; exact boundary TBD |
| Channel montage helpers | Electrode / parcel → atlas label mapping for consistent parcellation |
| Directed connectivity | Phase-slope index (PSI) or Granger causality; likely deferred to a separate package given the algorithmic distance from the current codebase |

------------------------------------------------------------------------

## Long-term

| Item | Notes |
|----|----|
| CRAN submission | Post full parity validation of all algorithms against Python dynfc |
| syncR integration | [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md) / [`stack_synchrony()`](https://dynr.circadia-lab.uk/reference/stack_synchrony.md) outputs → `syncR::sync()` participant database |
| Python parity vignette | End-to-end comparison dynR vs dynfc on the same dataset |

------------------------------------------------------------------------

## What stays out of scope

- **EEG preprocessing** (re-referencing, ICA, artifact rejection) →
  mrpheus
- **Brain state metrics** (fractional occupancy, dwell time, Markov) →
  stateR
- **Phase-amplitude coupling** — niche enough to be its own package
- **Autoregressive / Granger models** — different algorithmic family
