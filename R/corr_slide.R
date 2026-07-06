#' Sliding window correlation
#'
#' Compute functional connectivity matrices over sliding windows of a BOLD
#' timeseries (Hansen et al., 2015).
#'
#' The correlation matrices are computed by a compiled C++ backend
#' ([corr_slide_cpp()]) that calculates Pearson correlation directly on the
#' input matrix. No per-window submatrix allocation or transposition is
#' performed; the t-outer loop ordering exploits the column-major memory
#' layout of the input for sequential cache access.
#'
#' @param timeseries Numeric matrix \[N x Tmax\]. BOLD signal with N parcels
#'   as rows and Tmax timepoints as columns.
#' @param window Integer. Window size in timepoints.
#' @param step Integer. Step between window onsets. Defaults to `window`
#'   (non-overlapping windows).
#'
#' @return A list with:
#'   \item{corr_mats}{Array \[N, N, n_windows\] of Pearson correlation matrices.}
#'   \item{idx}{Integer vector of 1-indexed window onset positions.}
#'
#' @references
#' Hansen, E. C. A. et al. (2015). Functional connectivity dynamics: Modeling
#' the switching behavior of the resting state. *NeuroImage*, 105, 525-535.
#' \doi{10.1016/j.neuroimage.2014.11.001}
#'
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' res <- corr_slide(ts, window = 20)
#' dim(res$corr_mats)  # 10 x 10 x 10
corr_slide <- function(timeseries, window, step = NULL) {
  if (is.null(step)) step <- window
  corr_slide_cpp(timeseries, window, as.integer(step))
}
