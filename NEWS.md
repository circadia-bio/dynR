# NEWS.md

## dynR 0.1.2

### New functions

* `dyn_transitions()`: first-order Markov transition probabilities between
  brain states. Takes a long-format timepoint-level data frame, computes
  source → target transitions within each group (subject × session), and
  returns a nested tibble of per-transition probabilities ready for
  `grouped_perm_glmm()`. Ported from `clusters_markov()` in the
  `neonatal_dfc` analysis pipeline (França et al., *Nat Commun*).
  `remIntra = TRUE` strips self-transitions before normalising.

### Dependencies

* Added `dplyr`, `tidyr`, and `rlang` to `Imports` (required by
  `dyn_transitions()`).

---

## dynR 0.1.1

### Bug fixes

* `bandpass_filter()`: replaced `gsignal::filtfilt()` (zero initial conditions)
  with a scipy-compatible implementation that achieves bit-perfect parity with
  `scipy.signal.filtfilt` and `scipy.signal.sosfiltfilt` (max |diff| < 4e-12).
  The old implementation produced edge transients of up to ~0.24 signal units
  on real fMRI data due to zero initial state; this propagated into
  `hilbert_phases()` and `dyn_phase_lock()`, causing LEiDA eigenvectors to
  diverge from the Python `dynfc` reference pipeline.

  The fix introduces three internal helpers:
  - `.lfilter_zi()`: steady-state initial conditions via companion-matrix solve
    (equivalent to `scipy.signal.lfilter_zi`).
  - `.lfilter()`: Direct Form II Transposed IIR filter with explicit initial
    state (equivalent to `scipy.signal.lfilter`).
  - `.odd_ext()`: odd-reflection signal extension matching
    `scipy.signal._arraytools.odd_ext` exactly — right extension uses
    `x[N-2 : N-n-2 : -1]`, padlen = `3 × max(length(b), length(a))`.

  With these corrections, per-timepoint correlations between R and Python
  LEiDA eigenvectors are identically 1.0 across all sessions.

### Tests

* Added `tests/testthat/test-bandpass_filter.R` with six tests adapted from
  the `scipy.signal.filtfilt` documentation examples:
  output-length preservation, DC blocking, in-band sine-wave recovery,
  stopband attenuation (>40 dB), zero-phase property, and machine-precision
  agreement with `scipy.signal.filtfilt` on fMRI parameters (guarded by
  `skip_if_not_installed("reticulate")`).
* Added `reticulate` to `Suggests` for the scipy cross-validation test.

---

## dynR 0.1.0

* Initial scaffold.
* Ported core dynamic functional connectivity functions from the Python `dynfc` package.
* Implements: `bandpass_filter()`, `hilbert_phases()`, `corr_slide()`, `cofluct()`,
  `corr_corr()`, `dyn_phase_lock()`, `get_leida()`, `kuramoto()`, `shannon_entropy()`,
  `do_euclid()`.
