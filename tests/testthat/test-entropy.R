# Reference values reproduced from Python dynfc test suite (scipy.stats.entropy)
# series = [1, 3, 5, 2, 3, 5, 3, 2, 1, 3, 4, 5]
# getEntropy(series, base=None) == 1.5171063970610277  (natural log)
# getEntropy(series, base=2)    == 2.1887218755408675  (base-2)

test_that("shannon_entropy: natural log matches Python scipy reference", {
  series <- c(1, 3, 5, 2, 3, 5, 3, 2, 1, 3, 4, 5)
  expect_equal(
    shannon_entropy(series, base = exp(1)),
    1.5171063970610277,
    tolerance = 1e-10
  )
})

test_that("shannon_entropy: base-2 matches Python scipy reference", {
  series <- c(1, 3, 5, 2, 3, 5, 3, 2, 1, 3, 4, 5)
  expect_equal(
    shannon_entropy(series, base = 2),
    2.1887218755408675,
    tolerance = 1e-10
  )
})

test_that("shannon_entropy: uniform distribution gives maximum entropy", {
  # For n equally probable outcomes: H = log2(n)
  series <- 1:8  # 8 unique values, each count = 1
  expect_equal(shannon_entropy(series, base = 2), log2(8), tolerance = 1e-10)
})

test_that("shannon_entropy: constant series gives zero entropy", {
  series <- rep(3, 20)
  expect_equal(shannon_entropy(series), 0)
})

test_that("shannon_entropy: n_bits discretisation matches manual round-and-scale", {
  set.seed(1)
  x      <- runif(100)
  h_raw  <- shannon_entropy(round(x * 256), base = 2) # manual: scale then round
  h_bits <- shannon_entropy(x, base = 2, n_bits = 8)  # via n_bits
  expect_equal(h_raw, h_bits, tolerance = 1e-10)
})
