# Changelog

## dynR 0.1.5

### New functions

- [`leida_pipeline()`](https://dynr.circadia-lab.uk/reference/leida_pipeline.md)
  /
  [`sw_pipeline()`](https://dynr.circadia-lab.uk/reference/sw_pipeline.md):
  high-level wrappers running the full phase-based or correlation-based
  pipeline in one call. Both accept an optional `filter` argument to
  skip bandpass filtering on pre-filtered data. Return S3 objects
  (`dynR_leida` / `dynR_sw`) with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

- [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md): FC
  matrix heatmap using the dynR diverging palette (deep indigo -\>
  periwinkle -\> brick red).

- [`plot_synchrony()`](https://dynr.circadia-lab.uk/reference/plot_synchrony.md):
  Kuramoto R(t) time series with mean reference line.

- [`plot_state_sequence()`](https://dynr.circadia-lab.uk/reference/plot_state_sequence.md):
  brain state tile plot with dynR state palette.

- [`plot.dynR_leida()`](https://dynr.circadia-lab.uk/reference/plot.dynR_leida.md)
  /
  [`plot.dynR_sw()`](https://dynr.circadia-lab.uk/reference/plot.dynR_sw.md):
  S3 plot methods for pipeline outputs; dispatch to
  [`plot_synchrony()`](https://dynr.circadia-lab.uk/reference/plot_synchrony.md)
  / [`plot_fc()`](https://dynr.circadia-lab.uk/reference/plot_fc.md) or
  the RSS plot via `type = c("synchrony"|"fc")` and
  `type = c("rss"|"fc")`.

- [`batch_leida()`](https://dynr.circadia-lab.uk/reference/batch_leida.md)
  / [`batch_sw()`](https://dynr.circadia-lab.uk/reference/batch_sw.md):
  apply the respective pipeline to multiple subjects at once. Accepts
  either a named list of \[N x Tmax\] matrices or a 3-D array \[N x Tmax
  x subjects\].

- [`stack_leida()`](https://dynr.circadia-lab.uk/reference/stack_leida.md):
  combine LEiDA eigenvector matrices from a batch result into a single
  data frame or matrix for cross-subject K-means clustering.

- [`stack_synchrony()`](https://dynr.circadia-lab.uk/reference/stack_synchrony.md):
  combine synchrony time series into a tidy long-format data frame with
  `subject`, `timepoint`, `synchrony`, and `metastability` columns.

### Dependencies

- `ggplot2` moved from `Suggests` to `Imports` (required by plot
  functions).

------------------------------------------------------------------------

## dynR 0.1.4

### Performance

Four additional compiled backends and one vectorisation, completing the
C++ acceleration of all main compute paths:

- [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md):
  C++ / LAPACK backend
  ([`get_leida_cpp()`](https://dynr.circadia-lab.uk/reference/get_leida_cpp.md)).
  Calls LAPACK `dsyev` directly via `R_ext/Lapack.h`; workspace
  allocated once and reused across all timepoints. `dsyev` returns
  eigenvalues in ascending order so the leading eigenvector is the last
  column (opposite of R’s
  [`eigen()`](https://rdrr.io/r/base/eigen.html)). Parity vs
  [`eigen()`](https://rdrr.io/r/base/eigen.html): max diff \< 2.5e-15.

- [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md):
  C++ backend
  ([`kuramoto_sync_cpp()`](https://dynr.circadia-lab.uk/reference/kuramoto_sync_cpp.md)).
  Replaces `vapply` + R complex `exp(1i*x)` with direct `cos`/`sin`
  accumulation and `sqrt` — no complex-number allocation. Parity vs R
  `vapply` reference: \< 1e-14.

- [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md):
  vectorised with [`mvfft()`](https://rdrr.io/r/stats/fft.html).
  Replaces N per-parcel [`fft()`](https://rdrr.io/r/stats/fft.html)
  calls in an R loop with two
  [`mvfft()`](https://rdrr.io/r/stats/fft.html) calls on the full
  matrix. The Hilbert multiplier is broadcast across parcel columns in
  one step. (`fft_factor`/`fft_work` are not in R’s public API on
  macOS/Accelerate; `mvfft` achieves equivalent throughput without a
  bundled FFT.) Parity vs `.hilbert_r()` loop: \< 1e-14.

- [`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md):
  C++ backend
  ([`corr_slide_cpp()`](https://dynr.circadia-lab.uk/reference/corr_slide_cpp.md)).
  Computes Pearson correlation directly on the input matrix without
  per-window submatrix allocation or transposition. Cross-products
  accumulated in the upper triangle with t-outer loop ordering
  (column-major sequential access). Tested for non-overlapping and
  overlapping (step = 10, 5) windows. Parity vs
  [`cor()`](https://rdrr.io/r/stats/cor.html): \< 1e-10; removes
  [`stats::cor`](https://rdrr.io/r/stats/cor.html) from `Imports`.

- [`do_euclid()`](https://dynr.circadia-lab.uk/reference/do_euclid.md):
  vectorised with [`diff()`](https://rdrr.io/r/base/diff.html) +
  [`rowSums()`](https://rdrr.io/r/base/colSums.html). Loop over rows
  replaced by `c(0, sqrt(rowSums(diff(x)^2)))`.

### Tests

- Added parity tests for all new backends: `get_leida_cpp`,
  `kuramoto_sync_cpp`, `hilbert_phases` (mvfft path), `corr_slide_cpp`,
  and `dyn_phase_lock_cpp` (bit-perfect, max diff == 0).
- Added `tests/testthat/test-leida.R` covering dimensions, sign
  convention, and LAPACK/eigen parity.

------------------------------------------------------------------------

## dynR 0.1.3

### Performance

- [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md)
  now delegates the phase-locking loop to a compiled C++ backend
  ([`dyn_phase_lock_cpp()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock_cpp.md)).
  The inner `outer(cos(a-b))` R loop is replaced by a triple nested C++
  loop that exploits symmetry (`cos(a-b) = cos(b-a)`): only the upper
  triangle is evaluated and mirrored, halving trigonometric operations.
  Adds `Rcpp` as a dependency.

### Dependencies

- Removed `gsignal` from `Imports`. Both remaining uses have been ported
  to base-R equivalents:
  - `gsignal::butter()` replaced by `.butter_bandpass()` — a direct port
    of `scipy.signal.butter()` (analog Butterworth LP prototype -\>
    LP-to-BP transformation -\> bilinear transform -\> ZPK-to-TF).
    Coefficients validated against the `[k, 0, -2k, 0, k]` bandpass
    structure and existing scipy parity tests.
  - `gsignal::hilbert()` replaced by `.hilbert_r()` — FFT-based analytic
    signal using base-R [`fft()`](https://rdrr.io/r/stats/fft.html),
    equivalent to `scipy.signal.hilbert()`. Phase range confirmed
    `[-pi, pi]` on synthetic cosine.

------------------------------------------------------------------------

## dynR 0.1.2

### New functions

- [`dyn_transitions()`](https://dynr.circadia-lab.uk/reference/dyn_transitions.md):
  first-order Markov transition probabilities between brain states.
  Takes a long-format timepoint-level data frame, computes source →
  target transitions within each group (subject × session), and returns
  a nested tibble of per-transition probabilities ready for
  `grouped_perm_glmm()`. Ported from `clusters_markov()` in the
  `neonatal_dfc` analysis pipeline (França et al., *Nat Commun*).
  `remIntra = TRUE` strips self-transitions before normalising.

### Dependencies

- Added `dplyr`, `tidyr`, and `rlang` to `Imports` (required by
  [`dyn_transitions()`](https://dynr.circadia-lab.uk/reference/dyn_transitions.md)).

------------------------------------------------------------------------

## dynR 0.1.1

### Bug fixes

- [`bandpass_filter()`](https://dynr.circadia-lab.uk/reference/bandpass_filter.md):
  replaced `gsignal::filtfilt()` (zero initial conditions) with a
  scipy-compatible implementation that achieves bit-perfect parity with
  `scipy.signal.filtfilt` and `scipy.signal.sosfiltfilt` (max \|diff\|
  \< 4e-12). The old implementation produced edge transients of up to
  ~0.24 signal units on real fMRI data due to zero initial state; this
  propagated into
  [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md)
  and
  [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md),
  causing LEiDA eigenvectors to diverge from the Python `dynfc`
  reference pipeline.

  The fix introduces three internal helpers:

  - `.lfilter_zi()`: steady-state initial conditions via
    companion-matrix solve (equivalent to `scipy.signal.lfilter_zi`).
  - `.lfilter()`: Direct Form II Transposed IIR filter with explicit
    initial state (equivalent to `scipy.signal.lfilter`).
  - `.odd_ext()`: odd-reflection signal extension matching
    `scipy.signal._arraytools.odd_ext` exactly — right extension uses
    `x[N-2 : N-n-2 : -1]`, padlen = `3 × max(length(b), length(a))`.

  With these corrections, per-timepoint correlations between R and
  Python LEiDA eigenvectors are identically 1.0 across all sessions.

### Tests

- Added `tests/testthat/test-bandpass_filter.R` with six tests adapted
  from the `scipy.signal.filtfilt` documentation examples: output-length
  preservation, DC blocking, in-band sine-wave recovery, stopband
  attenuation (\>40 dB), zero-phase property, and machine-precision
  agreement with `scipy.signal.filtfilt` on fMRI parameters (guarded by
  `skip_if_not_installed("reticulate")`).
- Added `reticulate` to `Suggests` for the scipy cross-validation test.

------------------------------------------------------------------------

## dynR 0.1.0

- Initial scaffold.
- Ported core dynamic functional connectivity functions from the Python
  `dynfc` package.
- Implements:
  [`bandpass_filter()`](https://dynr.circadia-lab.uk/reference/bandpass_filter.md),
  [`hilbert_phases()`](https://dynr.circadia-lab.uk/reference/hilbert_phases.md),
  [`corr_slide()`](https://dynr.circadia-lab.uk/reference/corr_slide.md),
  [`cofluct()`](https://dynr.circadia-lab.uk/reference/cofluct.md),
  [`corr_corr()`](https://dynr.circadia-lab.uk/reference/corr_corr.md),
  [`dyn_phase_lock()`](https://dynr.circadia-lab.uk/reference/dyn_phase_lock.md),
  [`get_leida()`](https://dynr.circadia-lab.uk/reference/get_leida.md),
  [`kuramoto()`](https://dynr.circadia-lab.uk/reference/kuramoto.md),
  [`shannon_entropy()`](https://dynr.circadia-lab.uk/reference/shannon_entropy.md),
  [`do_euclid()`](https://dynr.circadia-lab.uk/reference/do_euclid.md).
