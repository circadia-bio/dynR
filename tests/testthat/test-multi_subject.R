# tests/testthat/test-multi_subject.R

make_ts_list <- function(n_sub = 2, N = 10, Tmax = 100, seed = 1) {
  set.seed(seed)
  nms <- paste0("sub", formatC(seq_len(n_sub), width = 2, flag = "0"))
  stats::setNames(
    lapply(seq_len(n_sub), function(i) matrix(rnorm(N * Tmax), nrow = N)),
    nms
  )
}

# ── batch_leida ───────────────────────────────────────────────────────────────

test_that("batch_leida: accepts a named list", {
  ts_list <- make_ts_list(2)
  res     <- batch_leida(ts_list, filter = FALSE)

  expect_type(res, "list")
  expect_length(res, 2L)
  expect_named(res, names(ts_list))
  expect_s3_class(res[[1L]], "dynR_leida")
})

test_that("batch_leida: accepts a 3-D array", {
  set.seed(2)
  arr <- array(rnorm(10 * 100 * 3), dim = c(10, 100, 3),
               dimnames = list(NULL, NULL, c("sub01", "sub02", "sub03")))
  res <- batch_leida(arr, filter = FALSE)

  expect_length(res, 3L)
  expect_named(res, c("sub01", "sub02", "sub03"))
  expect_s3_class(res[["sub01"]], "dynR_leida")
})

test_that("batch_leida: unnamed 3-D array gets synthetic subject names", {
  set.seed(3)
  arr <- array(rnorm(8 * 80 * 2), dim = c(8, 80, 2))
  res <- batch_leida(arr, filter = FALSE)

  expect_named(res, c("sub1", "sub2"))
})

test_that("batch_leida: passes extra arguments to leida_pipeline", {
  ts_list <- make_ts_list(2)
  res     <- batch_leida(ts_list, filter = FALSE)

  expect_equal(res[[1L]]$N,    nrow(ts_list[[1L]]))
  expect_equal(res[[1L]]$Tmax, ncol(ts_list[[1L]]))
})

# ── batch_sw ──────────────────────────────────────────────────────────────────

test_that("batch_sw: accepts a named list", {
  ts_list <- make_ts_list(2)
  res     <- batch_sw(ts_list, window = 20, filter = FALSE)

  expect_length(res, 2L)
  expect_s3_class(res[[1L]], "dynR_sw")
  expect_equal(res[[1L]]$window, 20L)
})

test_that("batch_sw: accepts a 3-D array", {
  set.seed(4)
  arr <- array(rnorm(10 * 100 * 2), dim = c(10, 100, 2),
               dimnames = list(NULL, NULL, c("A", "B")))
  res <- batch_sw(arr, window = 25, filter = FALSE)

  expect_named(res, c("A", "B"))
  expect_s3_class(res[["A"]], "dynR_sw")
})

# ── stack_leida ───────────────────────────────────────────────────────────────

test_that("stack_leida: returns data frame with subject column by default", {
  ts_list <- make_ts_list(3)
  batch   <- batch_leida(ts_list, filter = FALSE)
  stacked <- stack_leida(batch)

  expect_s3_class(stacked, "data.frame")
  expect_true("subject" %in% names(stacked))
  expect_equal(nrow(stacked), sum(sapply(batch, function(r) nrow(r$leida))))
})

test_that("stack_leida: add_subject_id = FALSE returns plain matrix", {
  ts_list <- make_ts_list(2)
  batch   <- batch_leida(ts_list, filter = FALSE)
  stacked <- stack_leida(batch, add_subject_id = FALSE)

  expect_true(is.matrix(stacked))
  expect_equal(ncol(stacked), nrow(ts_list[[1L]]))
})

test_that("stack_leida: stacked row count equals sum of per-subject frames", {
  ts_list <- make_ts_list(3)
  batch   <- batch_leida(ts_list, filter = FALSE)
  stacked <- stack_leida(batch)

  total_rows <- sum(sapply(batch, function(r) nrow(r$leida)))
  expect_equal(nrow(stacked), total_rows)
})

# ── stack_synchrony ───────────────────────────────────────────────────────────

test_that("stack_synchrony: returns tidy data frame with expected columns", {
  ts_list <- make_ts_list(2)
  batch   <- batch_leida(ts_list, filter = FALSE)
  df      <- stack_synchrony(batch)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("subject", "timepoint", "synchrony", "metastability")
                  %in% names(df)))
})

test_that("stack_synchrony: row count equals total timepoints across subjects", {
  ts_list <- make_ts_list(3)
  batch   <- batch_leida(ts_list, filter = FALSE)
  df      <- stack_synchrony(batch)

  total <- sum(sapply(batch, function(r) length(r$synchrony)))
  expect_equal(nrow(df), total)
})

test_that("stack_synchrony: metastability is constant within subject", {
  ts_list <- make_ts_list(2)
  batch   <- batch_leida(ts_list, filter = FALSE)
  df      <- stack_synchrony(batch)

  for (sub in unique(df$subject)) {
    m <- df$metastability[df$subject == sub]
    expect_true(length(unique(m)) == 1L)
  }
})
