# NEWS.md

## dynR 0.1.5

### New functions

* `leida_pipeline()` / `sw_pipeline()`: high-level wrappers running the full
  phase-based or correlation-based pipeline in one call. Both accept an
  optional `filter` argument to skip bandpass filtering on pre-filtered data.
  Return S3 objects (`dynR_leida` / `dynR_sw`) with `print()` and `plot()`
  methods.

* `plot_fc()`: FC matrix heatmap using the circadia diverging palette
  (deep blue -> antique white -> coral).

* `plot_synchrony()`: Kuramoto R(t) time series with mean reference line.

* `plot_state_sequence()`: brain state tile plot with circadia state colours.

* `plot.dynR_leida()` / `plot.dynR_sw()`: S3 plot methods for pipeline
  outputs; dispatch to `plot_synchrony()` / `plot_fc()` or the RSS plot
  via `type = c("synchrony"|"fc")` and `type = c("rss"|"fc")`.

* `batch_leida()` / `batch_sw()`: apply the respective pipeline to multiple
  subjects at once. Accepts either a named list of [N x Tmax] matrices or a
  3-D array [N x Tmax x subjects].

* `stack_leida()`: combine LEiDA eigenvector matrices from a batch result into
  a single data frame or matrix for cross-subject K-means clustering.

* `stack_synchrony()`: combine synchrony time series into a tidy long-format
  data frame with `subject`, `timepoint`, `synchrony`, and `metastability`
  columns.

### Dependencies

* `ggplot2` moved from `Suggests` to `Imports` (required by plot functions).

---

## dynR 0.1.4

### Performance

Four additional compiled backends and one vectorisation, completing the C++
acceleration of all main compute paths:

* `get_leida()`: C++ / LAPACK backend (`get_leida_cpp()`). Calls LAPACK
  `dsyev` directly via `R_ext/Lapack.h`; workspace allocated once and reused
  across all timepoints. `dsyev` returns eigenvalues in ascending order so the
  leading eigenvector is the last column (opposite of R's `eigen()`). Parity
  vs `eigen()`: max diff < 2.5e-15.

* `kuramoto()`: C++ backend (`kuramoto_sync_cpp()`). Replaces `vapply` +
  R complex `exp(1i*x)` with direct `cos`/`sin` accumulation and `sqrt` —
  no complex-number allocation. Parity vs R `vapply` reference: < 1e-14.

* `hilbert_phases()`: vectorised with `mvfft()`. Replaces N per-parcel
  `fft()` calls in an R loop with two `mvfft()` calls on the full matrix.
  The Hilbert multiplier is broadcast across parcel columns in one step.
  (`fft_factor`/`fft_work` are not in R's public API on macOS/Accelerate;
  `mvfft` achieves equivalent throughput without a bundled FFT.)
  Parity vs `.hilbert_r()` loop: < 1e-14.

* `corr_slide()`: C++ backend (`corr_slide_cpp()`). Computes Pearson
  correlation directly on the input matrix without per-window submatrix
  allocation or transposition. Cross-products accumulated in the upper
  triangle with t-outer loop ordering (column-major sequential access).
  Tested for non-overlapping and overlapping (step = 10, 5) windows.
  Parity vs `cor()`: < 1e-10; removes `stats::cor` from `Imports`.

* `do_euclid()`: vectorised with `diff()` + `rowSums()`. Loop over rows
  replaced by `c(0, sqrt(rowSums(diff(x)^2)))`.

### Tests

* Added parity tests for all new backends: `get_leida_cpp`, `kuramoto_sync_cpp`,
  `hilbert_phases` (mvfft path), `corr_slide_cpp`, and `dyn_phase_lock_cpp`
  (bit-perfect, max diff == 0).
* Added `tests/testthat/test-leida.R` covering dimensions, sign convention,
  and LAPACK/eigen parity.

---

## dynR 0.1.3

### Performance

* `dyn_phase_lock()` now delegates the phase-locking loop to a compiled C++
  backend (`dyn_phase_lock_cpp()`). The inner `outer(cos(a-b))` R loop is
  replaced by a triple nested C++ loop that exploits symmetry
  (`cos(a-b) = cos(b-a)`): only the upper triangle is evaluated and mirrored,
  halving trigonometric operations. Adds `Rcpp` as a dependency.

### Dependencies

* Removed `gsignal` from `Imports`. Both remaining uses have been ported to
  base-R equivalents:
  - `gsignal::butter()` replaced by `.butter_bandpass()` — a direct port of
    `scipy.signal.butter()` (analog Butterworth LP prototype -> LP-to-BP
    transformation -> bilinear transform -> ZPK-to-TF). Coefficients validated
    against the `[k, 0, -2k, 0, k]` bandpass structure and existing scipy
    parity tests.
  - `gsignal::hilbert()` replaced by `.hilbert_r()` — FFT-based analytic
    signal using base-R `fft()`, equivalent to `scipy.signal.hilbert()`.
    Phase range confirmed `[-pi, pi]` on synthetic cosine.

---

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
