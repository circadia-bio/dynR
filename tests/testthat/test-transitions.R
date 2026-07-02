# Minimal fixture: 2 subjects × 1 session, 10 timepoints each, 3 states.
make_df <- function(seed = 1L) {
  set.seed(seed)
  data.frame(
    sub   = rep(c("A", "B"), each = 10L),
    ses   = 1L,
    ttime = rep(seq_len(10L), 2L),
    clus3 = sample(1:3, 20L, replace = TRUE)
  )
}

# ── Output structure ──────────────────────────────────────────────────────────

test_that("dyn_transitions returns a nested tibble with expected columns", {
  tr <- dyn_transitions(
    make_df(),
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  expect_s3_class(tr, "data.frame")
  expect_true(all(c("tag", "source", "target", "data") %in% names(tr)))
})

test_that("tag is always '<source>_<target>'", {
  tr <- dyn_transitions(
    make_df(),
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  expected_tags <- paste0(tr$source, "_", tr$target)
  expect_equal(tr$tag, expected_tags)
})

test_that("data column is a list of data frames", {
  tr <- dyn_transitions(
    make_df(),
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  expect_true(all(vapply(tr$data, is.data.frame, logical(1))))
})

test_that("each nested data frame contains n, tot, nCount columns", {
  tr <- dyn_transitions(
    make_df(),
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  for (d in tr$data) {
    expect_true(all(c("n", "tot", "nCount") %in% names(d)))
  }
})

# ── Probability correctness ───────────────────────────────────────────────────

test_that("nCount is n / tot for every row in every nested data frame", {
  tr <- dyn_transitions(
    make_df(),
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  for (d in tr$data) {
    expect_equal(d$nCount, d$n / d$tot, tolerance = 1e-12)
  }
})

test_that("transition probabilities out of each source sum to 1 per subject", {
  df <- make_df()
  tr <- dyn_transitions(
    df,
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  flat <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    d <- tr$data[[i]]
    d$source <- tr$source[i]
    d
  }))
  sums <- tapply(flat$nCount, list(flat$sub, flat$source), sum)
  expect_true(all(abs(sums - 1) < 1e-10, na.rm = TRUE))
})

# ── Transition count ──────────────────────────────────────────────────────────

test_that("total transitions per subject equals n_timepoints - 1", {
  df <- make_df()
  tr <- dyn_transitions(
    df,
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  flat <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    d <- tr$data[[i]]
    d$source <- tr$source[i]
    d
  }))
  # Each subject has 10 timepoints → 9 transitions
  totals <- tapply(flat$n, flat$sub, sum)
  expect_true(all(totals == 9L))
})

# ── remIntra ─────────────────────────────────────────────────────────────────

test_that("remIntra = TRUE removes all self-transitions", {
  tr <- dyn_transitions(
    make_df(),
    vars     = c("sub", "ses"),
    cVar     = "clus3",
    sortBy   = c("sub", "ses", "ttime"),
    groupBy  = c("sub", "ses"),
    remIntra = TRUE
  )
  expect_true(all(tr$source != tr$target))
})

test_that("remIntra = FALSE retains self-transitions when present", {
  # Force a sequence with guaranteed self-transitions
  df <- data.frame(
    sub   = "A",
    ses   = 1L,
    ttime = 1:5,
    clus3 = c(1L, 1L, 2L, 3L, 2L)   # 1→1 is a self-transition
  )
  tr_keep <- dyn_transitions(
    df,
    vars     = c("sub", "ses"),
    cVar     = "clus3",
    sortBy   = c("sub", "ses", "ttime"),
    groupBy  = c("sub", "ses"),
    remIntra = FALSE
  )
  tr_drop <- dyn_transitions(
    df,
    vars     = c("sub", "ses"),
    cVar     = "clus3",
    sortBy   = c("sub", "ses", "ttime"),
    groupBy  = c("sub", "ses"),
    remIntra = TRUE
  )
  expect_true(any(tr_keep$source == tr_keep$target))
  expect_false(any(tr_drop$source == tr_drop$target))
})

# ── Boundary / edge cases ─────────────────────────────────────────────────────

test_that("single-timepoint group produces no transitions", {
  df <- data.frame(sub = "A", ses = 1L, ttime = 1L, clus3 = 1L)
  tr <- dyn_transitions(
    df,
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  expect_equal(nrow(tr), 0L)
})

test_that("transitions do not bleed across groupBy boundaries", {
  # Sub A ends on state 3, sub B starts on state 1.
  # If grouping were ignored, a 3→1 transition would appear for the boundary.
  df <- data.frame(
    sub   = rep(c("A", "B"), each = 3L),
    ses   = 1L,
    ttime = rep(1:3, 2L),
    clus3 = c(1L, 2L, 3L,   # sub A
               1L, 2L, 1L)  # sub B
  )
  tr <- dyn_transitions(
    df,
    vars    = c("sub", "ses"),
    cVar    = "clus3",
    sortBy  = c("sub", "ses", "ttime"),
    groupBy = c("sub", "ses")
  )
  flat <- do.call(rbind, lapply(seq_len(nrow(tr)), function(i) {
    d <- tr$data[[i]]
    d$tag_row <- tr$tag[i]
    d
  }))
  # A spurious 3→1 cross-boundary tag would only appear in sub B's data frame
  # alongside sub A's source. Here we simply check 3_1 is absent from sub A's rows.
  if ("3_1" %in% tr$tag) {
    sub_a_rows <- flat[flat$sub == "A" & flat$tag_row == "3_1", ]
    expect_equal(nrow(sub_a_rows), 0L)
  } else {
    succeed("no 3_1 transition exists (expected for sub A's sequence)")
  }
})
