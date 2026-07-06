test_that("kuramoto: synchrony length is Tmax - 20", {
  set.seed(1)
  N    <- 10L
  Tmax <- 100L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)
  ph   <- hilbert_phases(ts)
  res  <- kuramoto(ph)

  expect_equal(length(res$synchrony), Tmax - 20L)
})

test_that("kuramoto: synchrony values are in [0, 1]", {
  set.seed(2)
  ts  <- matrix(rnorm(8 * 120), nrow = 8)
  ph  <- hilbert_phases(ts)
  res <- kuramoto(ph)

  expect_true(all(res$synchrony >= 0 - 1e-10))
  expect_true(all(res$synchrony <= 1 + 1e-10))
})

test_that("kuramoto: metastability is a non-negative scalar", {
  set.seed(3)
  ts  <- matrix(rnorm(6 * 80), nrow = 6)
  ph  <- hilbert_phases(ts)
  res <- kuramoto(ph)

  expect_length(res$metastability, 1L)
  expect_gte(res$metastability, 0)
})

test_that("kuramoto: entropy is a finite non-negative scalar", {
  set.seed(4)
  ts  <- matrix(rnorm(5 * 60), nrow = 5)
  ph  <- hilbert_phases(ts)
  res <- kuramoto(ph)

  expect_length(res$entropy, 1L)
  expect_true(is.finite(res$entropy))
  expect_gte(res$entropy, 0)
})

test_that("kuramoto: perfectly synchronised signal gives synchrony near 1", {
  # All parcels with the same sinusoidal signal -> Kuramoto R ~= 1
  t_seq  <- seq(0, 10 * pi, length.out = 300)
  signal <- sin(t_seq)
  ts     <- matrix(rep(signal, 10), nrow = 10, byrow = TRUE)
  ph     <- hilbert_phases(ts)
  res    <- kuramoto(ph)

  expect_true(mean(res$synchrony) > 0.99)
})

test_that("kuramoto_sync_cpp: matches R vapply reference to machine precision", {
  # Confirms the C++ cos/sin accumulation agrees with R's complex-number
  # vapply implementation.
  set.seed(7)
  N    <- 20L
  Tmax <- 100L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)
  ph   <- hilbert_phases(ts)

  sync_cpp <- kuramoto(ph)$synchrony

  # R reference
  T_idx    <- seq(11L, Tmax - 10L)
  sync_ref <- vapply(T_idx, function(t) {
    Mod(sum(exp(1i * ph[, t])) / N)
  }, numeric(1))

  expect_equal(sync_cpp, sync_ref, tolerance = 1e-14)
})
