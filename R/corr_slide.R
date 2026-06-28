#' Sliding window correlation
#'
#' Compute functional connectivity matrices over sliding windows of a BOLD
#' timeseries (Hansen et al., 2015).
#'
#' @param timeseries Numeric matrix \[N × Tmax\]. BOLD signal with N parcels
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
#' the switching behavior of the resting state. *NeuroImage*, 105, 525–535.
#' \doi{10.1016/j.neuroimage.2014.11.001}
#'
#' @importFrom stats cor
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' res <- corr_slide(ts, window = 20)
#' dim(res$corr_mats)  # 10 x 10 x 10
corr_slide <- function(timeseries, window, step = NULL) {
  if (is.null(step)) step <- window
  N    <- nrow(timeseries)
  Tmax <- ncol(timeseries)
  idx  <- seq(1, Tmax, by = step)
  idx  <- idx[idx + window - 1 <= Tmax]
  n_windows <- length(idx)
  corr_mats <- array(0, dim = c(N, N, n_windows))
  for (w in seq_len(n_windows)) {
    seg <- timeseries[, idx[w]:(idx[w] + window - 1), drop = FALSE]
    corr_mats[, , w] <- cor(t(seg))
  }
  list(corr_mats = corr_mats, idx = idx)
}
