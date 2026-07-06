# tests/testthat/test-plot_dynR.R

# ── plot_fc ───────────────────────────────────────────────────────────────────

test_that("plot_fc: returns a ggplot object", {
  data(fc, package = "dynR", envir = environment())
  p <- plot_fc(fc[1:10, 1:10])
  expect_s3_class(p, "ggplot")
})

test_that("plot_fc: accepts custom title and limits", {
  m <- matrix(runif(25, -1, 1), 5, 5)
  p <- plot_fc(m, title = "Test FC", limits = c(-0.5, 0.5))
  expect_s3_class(p, "ggplot")
})

# ── plot_synchrony ────────────────────────────────────────────────────────────

test_that("plot_synchrony: returns a ggplot object", {
  set.seed(1)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  ph  <- hilbert_phases(ts)
  kop <- kuramoto(ph)
  p   <- plot_synchrony(kop$synchrony)
  expect_s3_class(p, "ggplot")
})

# ── plot_state_sequence ───────────────────────────────────────────────────────

test_that("plot_state_sequence: returns a ggplot object", {
  set.seed(2)
  states <- sample(1:4, 80, replace = TRUE)
  p <- plot_state_sequence(states)
  expect_s3_class(p, "ggplot")
})

test_that("plot_state_sequence: accepts factor input", {
  states <- factor(rep(c("A", "B", "C"), 10))
  p <- plot_state_sequence(states)
  expect_s3_class(p, "ggplot")
})

test_that("plot_state_sequence: accepts custom palette", {
  states <- sample(1:3, 30, replace = TRUE)
  p <- plot_state_sequence(states, palette = c("red", "blue", "green"))
  expect_s3_class(p, "ggplot")
})

# ── S3 plot methods ───────────────────────────────────────────────────────────

test_that("plot.dynR_leida: type = 'synchrony' returns ggplot", {
  set.seed(3)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- leida_pipeline(ts, filter = FALSE)
  p   <- plot(res, type = "synchrony")
  expect_s3_class(p, "ggplot")
})

test_that("plot.dynR_leida: type = 'fc' returns ggplot", {
  set.seed(4)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- leida_pipeline(ts, filter = FALSE)
  p   <- plot(res, type = "fc")
  expect_s3_class(p, "ggplot")
})

test_that("plot.dynR_sw: type = 'rss' returns ggplot", {
  set.seed(5)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- sw_pipeline(ts, window = 20, filter = FALSE)
  p   <- plot(res, type = "rss")
  expect_s3_class(p, "ggplot")
})

test_that("plot.dynR_sw: type = 'fc' returns ggplot", {
  set.seed(6)
  ts  <- matrix(rnorm(8 * 80), nrow = 8)
  res <- sw_pipeline(ts, window = 20, filter = FALSE)
  p   <- plot(res, type = "fc")
  expect_s3_class(p, "ggplot")
})
