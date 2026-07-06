#' Batch LEiDA pipeline across subjects
#'
#' Apply `leida_pipeline()` to multiple subjects in one call. Accepts either a
#' named list of \[N x Tmax\] matrices or a 3-D array \[N x Tmax x subjects\].
#'
#' @param data Named list of \[N x Tmax\] numeric matrices, **or** a 3-D
#'   numeric array with dimensions \[N, Tmax, subjects\]. For arrays, subject
#'   names are taken from `dimnames(data)[[3]]` if set, otherwise `"sub1"`,
#'   `"sub2"`, etc.
#' @param ... Additional arguments passed to `leida_pipeline()` (e.g. `flp`,
#'   `fhi`, `delt`, `filter`).
#'
#' @return Named list of `dynR_leida` objects, one per subject.
#' @seealso `stack_leida()`, `leida_pipeline()`
#' @export
#' @examples
#' set.seed(1)
#' ts_list <- list(
#'   sub01 = matrix(rnorm(10 * 200), nrow = 10),
#'   sub02 = matrix(rnorm(10 * 200), nrow = 10)
#' )
#' results <- batch_leida(ts_list, filter = FALSE)
#' names(results)
batch_leida <- function(data, ...) {
  data <- .coerce_to_ts_list(data)
  lapply(data, leida_pipeline, ...)
}

#' Batch sliding-window pipeline across subjects
#'
#' Apply `sw_pipeline()` to multiple subjects in one call. Accepts either a
#' named list of \[N x Tmax\] matrices or a 3-D array \[N x Tmax x subjects\].
#'
#' @param data Named list of \[N x Tmax\] matrices, or a 3-D array
#'   \[N, Tmax, subjects\]. See `batch_leida()`.
#' @param window Integer. Window size in timepoints passed to `sw_pipeline()`.
#' @param ... Additional arguments passed to `sw_pipeline()`.
#'
#' @return Named list of `dynR_sw` objects, one per subject.
#' @seealso `batch_leida()`, `sw_pipeline()`
#' @export
#' @examples
#' set.seed(2)
#' ts_list <- list(
#'   sub01 = matrix(rnorm(10 * 200), nrow = 10),
#'   sub02 = matrix(rnorm(10 * 200), nrow = 10)
#' )
#' results <- batch_sw(ts_list, window = 20, filter = FALSE)
#' names(results)
batch_sw <- function(data, window, ...) {
  data <- .coerce_to_ts_list(data)
  lapply(data, sw_pipeline, window = window, ...)
}

#' Stack LEiDA eigenvectors across subjects
#'
#' Combine the LEiDA eigenvector matrices from a list of `dynR_leida` results
#' into a single matrix (or data frame) ready for cross-subject K-means
#' clustering.
#'
#' @param batch_result Named list of `dynR_leida` objects, as returned by
#'   `batch_leida()`.
#' @param add_subject_id Logical. Prepend a `subject` column with the subject
#'   name? Default `TRUE`.
#'
#' @return If `add_subject_id = TRUE`: a data frame with columns `subject` and
#'   one column per parcel (named `V1`, `V2`, ...).
#'   If `FALSE`: a plain numeric matrix \[total_timepoints x N\].
#'
#' @seealso `batch_leida()`
#' @export
#' @examples
#' set.seed(3)
#' ts_list <- list(
#'   sub01 = matrix(rnorm(10 * 200), nrow = 10),
#'   sub02 = matrix(rnorm(10 * 200), nrow = 10)
#' )
#' results <- batch_leida(ts_list, filter = FALSE)
#' leida_all <- stack_leida(results)
#' dim(leida_all)
stack_leida <- function(batch_result, add_subject_id = TRUE) {
  nms  <- .subject_names(batch_result)
  mats <- mapply(function(res, nm) {
    m <- res$leida
    if (add_subject_id) {
      cbind(data.frame(subject = nm, stringsAsFactors = FALSE),
            as.data.frame(m))
    } else {
      m
    }
  }, batch_result, nms, SIMPLIFY = FALSE)
  do.call(rbind, mats)
}

#' Stack synchrony vectors across subjects
#'
#' Combine the Kuramoto synchrony time series from a list of `dynR_leida`
#' results into a tidy data frame with one row per timepoint per subject.
#'
#' @param batch_result Named list of `dynR_leida` objects.
#'
#' @return A data frame with columns:
#'   \item{subject}{Subject identifier.}
#'   \item{timepoint}{Timepoint index (within subject).}
#'   \item{synchrony}{Kuramoto R(t).}
#'   \item{metastability}{Subject-level metastability (constant within subject).}
#'
#' @seealso `batch_leida()`
#' @export
#' @examples
#' set.seed(4)
#' ts_list <- list(
#'   sub01 = matrix(rnorm(10 * 200), nrow = 10),
#'   sub02 = matrix(rnorm(10 * 200), nrow = 10)
#' )
#' results <- batch_leida(ts_list, filter = FALSE)
#' df_sync <- stack_synchrony(results)
#' head(df_sync)
stack_synchrony <- function(batch_result) {
  nms  <- .subject_names(batch_result)
  rows <- mapply(function(res, nm) {
    data.frame(
      subject          = nm,
      timepoint        = seq_along(res$synchrony),
      synchrony        = res$synchrony,
      metastability    = res$metastability,
      stringsAsFactors = FALSE
    )
  }, batch_result, nms, SIMPLIFY = FALSE)
  do.call(rbind, rows)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

# Coerce a 3-D array [N, Tmax, subjects] to a named list of [N x Tmax] matrices.
.coerce_to_ts_list <- function(data) {
  if (is.array(data) && length(dim(data)) == 3L) {
    n_sub <- dim(data)[3L]
    nms   <- if (!is.null(dimnames(data)[[3L]])) dimnames(data)[[3L]] else
               paste0("sub", seq_len(n_sub))
    data  <- stats::setNames(
      lapply(seq_len(n_sub), function(i) data[, , i]),
      nms
    )
  }
  if (!is.list(data))
    stop("'data' must be a named list of matrices or a 3-D array.", call. = FALSE)
  data
}

# Extract or synthesise subject names from a batch result list.
.subject_names <- function(batch_result) {
  nms <- names(batch_result)
  if (is.null(nms)) paste0("sub", seq_along(batch_result)) else nms
}
