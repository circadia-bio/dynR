test_that("do_euclid: first element is zero", {
  set.seed(1)
  x <- matrix(rnorm(50 * 3), nrow = 50)
  d <- do_euclid(x)

  expect_equal(d[1], 0)
})

test_that("do_euclid: output length equals nrow(input)", {
  set.seed(2)
  x <- matrix(rnorm(20 * 5), nrow = 20)
  d <- do_euclid(x)

  expect_length(d, nrow(x))
})

test_that("do_euclid: all distances are non-negative", {
  set.seed(3)
  x <- matrix(rnorm(30 * 4), nrow = 30)
  d <- do_euclid(x)

  expect_true(all(d >= 0))
})

test_that("do_euclid: known distance is correct", {
  # Two points: (0,0) and (3,4) → distance = 5
  x <- matrix(c(0, 0, 3, 4), nrow = 2, byrow = TRUE)
  d <- do_euclid(x)

  expect_equal(d[2], 5)
})

test_that("do_euclid: single-column matrix reduces to absolute differences", {
  x <- matrix(c(1, 4, 7, 3), ncol = 1)
  d <- do_euclid(x)

  expect_equal(d, c(0, 3, 3, 4))
})
