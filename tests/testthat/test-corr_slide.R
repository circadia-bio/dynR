test_that("corr_slide: single window equals full-timeseries correlation matrix", {
  data(ts, package = "dynR", envir = environment())
  data(fc, package = "dynR", envir = environment())

  res <- corr_slide(ts, window = ncol(ts))

  expect_equal(length(res$idx), 1L)
  expect_equal(res$corr_mats[, , 1], fc, tolerance = 1e-10)
})

test_that("corr_slide: two non-overlapping windows on doubled timeseries", {
  data(ts, package = "dynR", envir = environment())
  data(fc, package = "dynR", envir = environment())

  ts_double <- cbind(ts, ts)
  res <- corr_slide(ts_double, window = ncol(ts))

  expect_equal(dim(res$corr_mats)[3], 2L)
  diff <- sum(abs(res$corr_mats[, , 1] - fc)) +
          sum(abs(res$corr_mats[, , 2] - fc))
  expect_equal(round(diff), 0)
})

test_that("corr_slide: output dimensions are correct", {
  set.seed(1)
  ts <- matrix(rnorm(10 * 100), nrow = 10)
  res <- corr_slide(ts, window = 20, step = 10)

  expect_equal(dim(res$corr_mats)[1], 10L)
  expect_equal(dim(res$corr_mats)[2], 10L)
  expect_equal(dim(res$corr_mats)[3], length(res$idx))
})

test_that("corr_slide: diagonal of each window matrix is 1", {
  set.seed(2)
  ts <- matrix(rnorm(5 * 60), nrow = 5)
  res <- corr_slide(ts, window = 20)

  for (w in seq_len(dim(res$corr_mats)[3])) {
    expect_equal(diag(res$corr_mats[, , w]), rep(1, 5), tolerance = 1e-10)
  }
})

test_that("corr_slide: matrices are symmetric", {
  set.seed(3)
  ts <- matrix(rnorm(8 * 80), nrow = 8)
  res <- corr_slide(ts, window = 20)

  for (w in seq_len(dim(res$corr_mats)[3])) {
    expect_equal(res$corr_mats[, , w], t(res$corr_mats[, , w]), tolerance = 1e-12)
  }
})
