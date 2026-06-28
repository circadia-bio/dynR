test_that("cofluct: edge time series matches Python dynfc reference (ets10)", {
  ts10  <- load_fixture("ts10.rds")
  ets10 <- load_fixture("ets10.rds")

  res <- cofluct(ts10)

  # sum of absolute differences should round to 0
  diff <- sum(abs(res$edge_ts - ets10))
  expect_equal(round(diff), 0)
})

test_that("cofluct: rss length equals number of timepoints", {
  ts10 <- load_fixture("ts10.rds")
  res  <- cofluct(ts10)

  expect_equal(length(res$rss), ncol(ts10))
})

test_that("cofluct: number of edges is N*(N-1)/2 for k=1", {
  set.seed(1)
  N  <- 8L
  ts <- matrix(rnorm(N * 60), nrow = N)
  res <- cofluct(ts)

  expect_equal(nrow(res$edge_ts), N * (N - 1L) / 2L)
})

test_that("cofluct: rss is non-negative", {
  set.seed(2)
  ts  <- matrix(rnorm(6 * 50), nrow = 6)
  res <- cofluct(ts)

  expect_true(all(res$rss >= 0))
})

test_that("corr_corr: output is a square symmetric matrix", {
  set.seed(3)
  N    <- 6L
  Tmax <- 40L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)
  cc   <- corr_corr(ts)

  expect_equal(dim(cc), c(Tmax, Tmax))
  expect_equal(cc, t(cc), tolerance = 1e-12)
})

test_that("corr_corr: diagonal is 1", {
  set.seed(4)
  ts <- matrix(rnorm(5 * 30), nrow = 5)
  cc <- corr_corr(ts)

  expect_equal(diag(cc), rep(1, 30), tolerance = 1e-10)
})
