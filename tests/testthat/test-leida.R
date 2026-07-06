# tests/testthat/test-leida.R
#
# Parity tests for get_leida() and its C++ / LAPACK backend.

test_that("get_leida_cpp: matches R eigen() to machine precision", {
  # Confirms that LAPACK dsyev and R's eigen() agree on the leading
  # eigenvector to within floating-point rounding.
  #
  # Key differences that are accounted for:
  #   - dsyev returns eigenvalues in ascending order (leading = last column)
  #   - R's eigen() returns eigenvalues in descending order (leading = first column)
  #   - Both apply the same sign convention: row sum forced <= 0
  set.seed(7)
  N    <- 20L
  Tmax <- 100L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)
  ph   <- hilbert_phases(ts)
  sc   <- dyn_phase_lock_cpp(ph)

  leida_cpp <- get_leida(sc)

  # R reference (eigen, descending order)
  t_pts     <- dim(sc)[3]
  leida_ref <- matrix(0, nrow = t_pts, ncol = N)
  for (i in seq_len(t_pts)) {
    ev <- eigen(sc[, , i], symmetric = TRUE)
    v1 <- ev$vectors[, 1]
    if (sum(v1) > 0) v1 <- -v1
    leida_ref[i, ] <- v1
  }

  expect_lt(max(abs(leida_cpp - leida_ref)), 1e-14)
})

test_that("get_leida: sign convention (row sums <= 0)", {
  set.seed(8)
  ts  <- matrix(rnorm(15 * 80), nrow = 15)
  ph  <- hilbert_phases(ts)
  sc  <- dyn_phase_lock_cpp(ph)

  leida <- get_leida(sc)

  expect_true(all(rowSums(leida) <= 1e-10))
})

test_that("get_leida: output dimensions are correct", {
  set.seed(9)
  N    <- 12L
  Tmax <- 80L
  ts   <- matrix(rnorm(N * Tmax), nrow = N)
  ph   <- hilbert_phases(ts)
  sc   <- dyn_phase_lock_cpp(ph)

  leida <- get_leida(sc)

  expect_equal(dim(leida), c(Tmax - 20L, N))
})
