# Phase difference: cos(a - b)
# Mirrors the three phDiff tests in the Python dynfc test suite.

test_that("phase difference: cos(pi - pi) == 1", {
  expect_equal(cos(pi - pi), 1)
})

test_that("phase difference: cos(pi - pi/2) == 0", {
  expect_equal(round(cos(pi - pi / 2), 1), 0)
})

test_that("phase difference: cos(pi - 0) == -1", {
  expect_equal(round(cos(pi - 0), 1), -1)
})

test_that("dyn_phase_lock: sync_conn dimensions are correct", {
  set.seed(1)
  N     <- 10L
  Tmax  <- 100L
  ts    <- matrix(rnorm(N * Tmax), nrow = N)
  ph    <- hilbert_phases(ts)
  res   <- dyn_phase_lock(ph)

  expect_equal(dim(res$sync_conn), c(N, N, Tmax - 20L))
})

test_that("dyn_phase_lock: leida dimensions are correct", {
  set.seed(2)
  N     <- 10L
  Tmax  <- 100L
  ts    <- matrix(rnorm(N * Tmax), nrow = N)
  ph    <- hilbert_phases(ts)
  res   <- dyn_phase_lock(ph)

  expect_equal(dim(res$leida), c(Tmax - 20L, N))
})

test_that("dyn_phase_lock: each sync_conn slice is symmetric", {
  set.seed(3)
  ts  <- matrix(rnorm(6 * 80), nrow = 6)
  ph  <- hilbert_phases(ts)
  res <- dyn_phase_lock(ph)

  for (t in seq_len(dim(res$sync_conn)[3])) {
    expect_equal(res$sync_conn[, , t], t(res$sync_conn[, , t]), tolerance = 1e-12)
  }
})

test_that("dyn_phase_lock: diagonal of each sync_conn slice is 1", {
  set.seed(4)
  N   <- 5L
  ts  <- matrix(rnorm(N * 60), nrow = N)
  ph  <- hilbert_phases(ts)
  res <- dyn_phase_lock(ph)

  for (t in seq_len(dim(res$sync_conn)[3])) {
    expect_equal(diag(res$sync_conn[, , t]), rep(1, N), tolerance = 1e-12)
  }
})

test_that("get_leida: sign convention (row sums <= 0)", {
  set.seed(5)
  ts     <- matrix(rnorm(8 * 80), nrow = 8)
  ph     <- hilbert_phases(ts)
  res    <- dyn_phase_lock(ph)

  row_sums <- rowSums(res$leida)
  expect_true(all(row_sums <= 1e-10))
})
