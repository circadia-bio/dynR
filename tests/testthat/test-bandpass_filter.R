# tests/testthat/test-bandpass_filter.R
#
# Tests for bandpass_filter() and its internal helpers.
#
# The numerical accuracy tests are inspired by the examples in the
# scipy.signal.filtfilt documentation:
#   https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.filtfilt.html
#
# scipy docs example 1 (adapted for bandpass below):
#   t     = linspace(0, 1.0, 2001);  fs = 2000 Hz
#   xlow  = sin(2π × 5 × t)          ← low-freq component
#   xhigh = sin(2π × 250 × t)        ← high-freq component
#   x     = xlow + xhigh
#   b, a  = butter(8, 0.125)          ← 8th-order lowpass, 125 Hz cutoff
#   y     = filtfilt(b, a, x, padlen=150)
#   max(|y − xlow|)  ≈  9.1e-6
#
# The bandpass version replaces the lowpass with butter(8, [1, 100], 'bandpass')
# at the same sampling rate; the in-band component (5 Hz) should again be
# recovered with sub-1e-3 error in the middle of the signal.

# ── helpers ───────────────────────────────────────────────────────────────────

make_sine <- function(freq, n, fs) {
  t <- seq(0, (n - 1) / fs, by = 1 / fs)
  sin(2 * pi * freq * t)
}

# ── 1. Output length ──────────────────────────────────────────────────────────

test_that("bandpass_filter: output length equals input length", {
  set.seed(1)
  x <- rnorm(500)
  y <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 1, order = 2)
  expect_equal(length(y), 500L)
})

# ── 2. DC blocking ────────────────────────────────────────────────────────────

test_that("bandpass_filter: constant signal is filtered to near-zero (DC gain = 0)", {
  # A Butterworth bandpass has exactly zero DC gain, so any constant input
  # must produce an output of essentially zero.
  x <- rep(3.7, 1000)
  y <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 1, order = 2)
  expect_true(max(abs(y)) < 1e-8)
})

# ── 3. In-band recovery (adapted from scipy filtfilt docs example 1) ──────────

test_that("bandpass_filter: recovers in-band sine from a mixed signal (scipy docs analogy)", {
  # Analogous to scipy docs example 1, but using a bandpass instead of a
  # lowpass. The original test used flp=1 Hz at fs=2000 Hz, which requires
  # ~2560 samples for the 8th-order filter to settle — longer than the entire
  # 2001-sample signal. Valid parameters need the lower cutoff high enough
  # that 1/f_lp << n/fs (filter impulse response fits inside the signal).
  #
  #   x_in  = sin(2pi x 200 x t)   <- 200 Hz, inside 100-600 Hz passband
  #   x_out = sin(2pi x 900 x t)   <- 900 Hz, outside passband
  #   x     = x_in + x_out
  #   bandpass order 8, cutoffs 100-600 Hz, fs = 2000 Hz
  #
  # scipy sosfiltfilt achieves max|y - x_in| < 1e-5 for these parameters.
  # We allow 1e-3 tolerance.

  n  <- 2001L
  fs <- 2000         # Hz
  x_in  <- make_sine(200, n, fs)   # 200 Hz, in passband
  x_out <- make_sine(900, n, fs)   # 900 Hz, out of passband
  x     <- x_in + x_out

  y   <- bandpass_filter(x, flp = 100, fhi = 600, delt = 1 / fs, order = 8)
  mid <- 200L:1800L

  expect_lt(max(abs(y[mid] - x_in[mid])), 1e-3)
})

# ── 4. Stopband attenuation ───────────────────────────────────────────────────

test_that("bandpass_filter: out-of-band signal is strongly attenuated (>40 dB)", {
  # A pure sine at 200 Hz is well outside the 1–10 Hz passband and should be
  # reduced by at least two orders of magnitude in RMS amplitude.
  n  <- 5001L
  fs <- 1000   # Hz
  x  <- make_sine(200, n, fs)  # 200 Hz, stopband

  y <- bandpass_filter(x, flp = 1, fhi = 10, delt = 1 / fs, order = 2)

  rms_ratio <- sqrt(mean(y^2)) / sqrt(mean(x^2))
  expect_lt(rms_ratio, 0.01)  # > 40 dB attenuation
})

# ── 5. Zero-phase property (filtfilt) ─────────────────────────────────────────

test_that("bandpass_filter: no phase shift on in-band cosine (zero-phase filtfilt)", {
  # filtfilt applies the filter forward then backward, giving zero net phase
  # shift. A pure in-band cosine should emerge with the same phase as the
  # input; only the amplitude may differ due to the passband gain.
  n  <- 2001L
  fs <- 100     # Hz
  f0 <- 3       # Hz; in passband 1–10 Hz
  t  <- seq(0, (n - 1) / fs, by = 1 / fs)
  x  <- cos(2 * pi * f0 * t)

  y   <- bandpass_filter(x, flp = 1, fhi = 10, delt = 1 / fs, order = 2)
  mid <- 200L:1800L

  # Normalise amplitudes and check that the zero-lag correlation ≈ 1
  y_mid <- y[mid] / sqrt(mean(y[mid]^2))
  x_mid <- x[mid] / sqrt(mean(x[mid]^2))
  expect_lt(max(abs(y_mid - x_mid)), 0.01)
})

# ── 6. fMRI parameters match scipy.signal.filtfilt to machine precision ───────

test_that("bandpass_filter: fMRI params match scipy.signal.filtfilt to 1e-9", {
  skip_if_not_installed("reticulate")
  skip_if(
    !tryCatch({ reticulate::import("scipy.signal"); TRUE }, error = function(e) FALSE),
    "scipy Python module not available"
  )

  set.seed(42)
  x   <- as.double(rnorm(515) * 1000 + 5000)   # fMRI-scale signal
  TR  <- 0.933
  flp <- 0.02
  fhi <- 0.10

  r_out <- bandpass_filter(x, flp = flp, fhi = fhi, delt = TR, order = 2L)

  sp  <- reticulate::import("scipy.signal")
  np  <- reticulate::import("numpy")
  ba  <- sp$butter(2L, c(flp, fhi), btype = "bandpass", fs = 1 / TR, output = "ba")
  b_r <- as.double(reticulate::py_to_r(ba[[1]]))
  a_r <- as.double(reticulate::py_to_r(ba[[2]]))

  py_out <- as.double(
    sp$filtfilt(
      np$array(b_r),
      np$array(a_r),
      np$array(as.double(x))
    )
  )

  expect_equal(r_out, py_out, tolerance = 1e-9)
})
