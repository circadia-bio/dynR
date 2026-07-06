# tests/testthat/test-hilbert_phases.R

test_that("hilbert_phases: output dimensions match input", {
  set.seed(1)
  ts     <- matrix(rnorm(10 * 200), nrow = 10)
  phases <- hilbert_phases(ts)

  expect_equal(dim(phases), dim(ts))
})

test_that("hilbert_phases: phases are in [-pi, pi]", {
  set.seed(2)
  ts     <- matrix(rnorm(8 * 150), nrow = 8)
  phases <- hilbert_phases(ts)

  expect_true(all(phases >= -pi - 1e-10))
  expect_true(all(phases <=  pi + 1e-10))
})

test_that("hilbert_phases: demeaning does not change phase structure", {
  # Shifted and unshifted signals should yield the same phases (Hilbert is
  # applied to demeaned signal in both cases)
  set.seed(3)
  row1     <- rnorm(100)
  ts_orig  <- matrix(row1, nrow = 1)
  ts_shift <- matrix(row1 + 5, nrow = 1)

  expect_equal(hilbert_phases(ts_orig), hilbert_phases(ts_shift), tolerance = 1e-10)
})

test_that("hilbert_phases: mvfft vectorised path matches .hilbert_r loop", {
  # Confirms the vectorised mvfft path (two calls on the full matrix) gives
  # identical results to the original per-parcel .hilbert_r() loop.
  # Both use R's FFT engine; any difference would be a sign/index bug.
  set.seed(7)
  N    <- 20L
  Tmax <- 100L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)

  phases_new <- hilbert_phases(ts)

  phases_ref <- matrix(0, nrow = N, ncol = Tmax)
  for (i in seq_len(N)) {
    x_c             <- ts[i, ] - mean(ts[i, ])
    phases_ref[i, ] <- Arg(dynR:::.hilbert_r(x_c))
  }

  expect_equal(phases_new, phases_ref, tolerance = 1e-14)
})
