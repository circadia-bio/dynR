# tests/testthat/test-pipeline.R

# ── leida_pipeline ────────────────────────────────────────────────────────────

test_that("leida_pipeline: returns dynR_leida object with correct structure", {
  set.seed(1)
  N <- 10L; Tmax <- 100L
  ts  <- matrix(rnorm(N * Tmax), nrow = N)
  res <- leida_pipeline(ts, filter = FALSE)

  expect_s3_class(res, "dynR_leida")
  expect_named(res, c("leida", "sync_conn", "phases", "synchrony",
                       "metastability", "entropy", "N", "Tmax"))
})

test_that("leida_pipeline: output dimensions are correct", {
  set.seed(2)
  N <- 10L; Tmax <- 100L
  ts  <- matrix(rnorm(N * Tmax), nrow = N)
  res <- leida_pipeline(ts, filter = FALSE)

  expect_equal(dim(res$leida),     c(Tmax - 20L, N))
  expect_equal(dim(res$sync_conn), c(N, N, Tmax - 20L))
  expect_equal(dim(res$phases),    c(N, Tmax))
  expect_equal(length(res$synchrony), Tmax - 20L)
  expect_equal(res$N,    N)
  expect_equal(res$Tmax, Tmax)
})

test_that("leida_pipeline: metastability and entropy are finite non-negative", {
  set.seed(3)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- leida_pipeline(ts, filter = FALSE)

  expect_true(is.finite(res$metastability) && res$metastability >= 0)
  expect_true(is.finite(res$entropy)       && res$entropy       >= 0)
})

test_that("leida_pipeline: print method works without error", {
  set.seed(4)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- leida_pipeline(ts, filter = FALSE)

  expect_output(print(res), "dynR_leida")
  expect_output(print(res), "Parcels")
})

test_that("leida_pipeline: matches individual function outputs", {
  set.seed(5)
  N <- 10L; Tmax <- 80L
  ts  <- matrix(rnorm(N * Tmax), nrow = N)
  res <- leida_pipeline(ts, filter = FALSE)

  ph  <- hilbert_phases(ts)
  dpl <- dyn_phase_lock(ph)
  kop <- kuramoto(ph)

  expect_equal(res$leida,         dpl$leida)
  expect_equal(res$synchrony,     kop$synchrony)
  expect_equal(res$metastability, kop$metastability)
})

# ── sw_pipeline ───────────────────────────────────────────────────────────────

test_that("sw_pipeline: returns dynR_sw object with correct structure", {
  set.seed(6)
  ts  <- matrix(rnorm(10 * 100), nrow = 10)
  res <- sw_pipeline(ts, window = 20, filter = FALSE)

  expect_s3_class(res, "dynR_sw")
  expect_named(res, c("corr_mats", "idx", "edge_ts", "rss",
                       "window", "step", "N", "Tmax"))
})

test_that("sw_pipeline: output dimensions are correct", {
  set.seed(7)
  N <- 10L; Tmax <- 100L; window <- 20L
  ts  <- matrix(rnorm(N * Tmax), nrow = N)
  res <- sw_pipeline(ts, window = window, step = window, filter = FALSE)

  n_windows <- floor(Tmax / window)
  expect_equal(dim(res$corr_mats)[1L], N)
  expect_equal(dim(res$corr_mats)[2L], N)
  expect_equal(dim(res$corr_mats)[3L], n_windows)
  expect_equal(length(res$rss), Tmax)
  expect_equal(res$window, window)
  expect_equal(res$step,   window)
})

test_that("sw_pipeline: print method works without error", {
  set.seed(8)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- sw_pipeline(ts, window = 20, filter = FALSE)

  expect_output(print(res), "dynR_sw")
  expect_output(print(res), "N windows")
})

test_that("sw_pipeline: matches individual function outputs", {
  set.seed(9)
  N <- 10L; Tmax <- 80L
  ts  <- matrix(rnorm(N * Tmax), nrow = N)
  res <- sw_pipeline(ts, window = 20, step = 10, filter = FALSE)

  sw <- corr_slide(ts, window = 20, step = 10)
  ec <- cofluct(ts)

  expect_equal(res$corr_mats, sw$corr_mats)
  expect_equal(res$idx,       sw$idx)
  expect_equal(res$rss,       ec$rss)
})
