#' Edge-centric cofluctuation analysis
#'
#' Compute edge time series and root-sum-square (RSS) cofluctuations from
#' z-standardised BOLD signal. Implements the edge-centric FC framework of
#' Esfahlani et al. (2020) and Faskowitz et al. (2020).
#'
#' @param timeseries Numeric matrix \[N × Tmax\]. BOLD signal with N parcels
#'   as rows and Tmax timepoints as columns.
#' @param k Integer. Upper-triangle offset. `k = 1` (default) excludes the
#'   diagonal; `k = 0` includes it.
#'
#' @return A list with:
#'   \item{edge_ts}{Numeric matrix \[n_edges × Tmax\]. Edge time series, one
#'     row per unique parcel pair.}
#'   \item{rss}{Numeric vector \[Tmax\]. Root-sum-square cofluctuation at
#'     each timepoint.}
#'
#' @references
#' Esfahlani, F. Z. et al. (2020). High-amplitude cofluctuations in cortical
#' activity drive functional connectivity. *PNAS*, 117(45), 28393–28401.
#' \doi{10.1073/pnas.2005531117}
#'
#' Faskowitz, J. et al. (2020). Edge-centric functional network representations
#' of human cerebral cortex reveal overlapping system-level architecture.
#' *Nature Neuroscience*, 23(12), 1644–1654.
#' \doi{10.1038/s41593-020-00719-y}
#'
#' @importFrom stats cor
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' res <- cofluct(ts)
#' dim(res$edge_ts)  # n_edges x 200
cofluct <- function(timeseries, k = 1) {
  N    <- nrow(timeseries)
  ts_z <- t(scale(t(timeseries)))  # z-score each parcel (row), ddof = 1
  idx  <- which(upper.tri(matrix(0L, N, N), diag = (k == 0)), arr.ind = TRUE)
  if (k > 1) idx <- idx[idx[, 2] - idx[, 1] >= k, , drop = FALSE]
  edge_ts <- ts_z[idx[, 1], ] * ts_z[idx[, 2], ]
  rss     <- sqrt(colSums(edge_ts^2))
  list(edge_ts = edge_ts, rss = rss)
}

#' Correlation of correlations matrix
#'
#' Compute the correlation between edge time series across all timepoints,
#' producing a \[Tmax × Tmax\] "correlation of correlations" matrix
#' (Hansen et al., 2015).
#'
#' @param timeseries Numeric matrix \[N × Tmax\]. BOLD signal.
#' @param k Integer. Upper-triangle offset passed to [cofluct()]. Default 1.
#'
#' @return Numeric matrix \[Tmax × Tmax\].
#'
#' @references
#' Hansen, E. C. A. et al. (2015). Functional connectivity dynamics: Modeling
#' the switching behavior of the resting state. *NeuroImage*, 105, 525–535.
#' \doi{10.1016/j.neuroimage.2014.11.001}
#'
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 100), nrow = 10, ncol = 100)
#' cc <- corr_corr(ts)
#' dim(cc)  # 100 x 100
corr_corr <- function(timeseries, k = 1) {
  res <- cofluct(timeseries, k)
  cor(res$edge_ts)
}
