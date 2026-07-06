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

test_that("corr_slide_cpp: matches R cor() reference to machine precision", {
  # Confirms the C++ Pearson implementation agrees with R's cor() to
  # well within machine precision. Tests both non-overlapping and
  # overlapping (step < window) cases.
  set.seed(7)
  N    <- 20L
  Tmax <- 120L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)

  for (step in c(20L, 10L, 5L)) {
    res_cpp <- corr_slide(ts, window = 20L, step = step)

    # R reference
    idx_r <- seq(1L, Tmax, by = step)
    idx_r <- idx_r[idx_r + 20L - 1L <= Tmax]
    ref   <- array(0, dim = c(N, N, length(idx_r)))
    for (w in seq_along(idx_r)) {
      seg        <- ts[, idx_r[w]:(idx_r[w] + 19L), drop = FALSE]
      ref[, , w] <- cor(t(seg))
    }

    expect_equal(res_cpp$corr_mats, ref,          tolerance = 1e-10)
    expect_equal(res_cpp$idx,       idx_r)
  }
})
