#' Hilbert transform phase extraction
#'
#' Compute the instantaneous phase time series for all parcels/voxels via
#' the analytic signal (Hilbert transform). Each parcel timeseries is
#' demeaned before transformation.
#'
#' @param timeseries Numeric matrix \[N × Tmax\]. BOLD signal with N parcels
#'   as rows and Tmax timepoints as columns.
#'
#' @return Numeric matrix \[N × Tmax\]. Instantaneous phases in radians.
#'
#' @references
#' Cabral, J. et al. (2017). Cognitive performance in healthy older adults
#' relates to spontaneous switching between states of functional connectivity
#' during rest. *Scientific Reports*, 7(1), 5135.
#' \doi{10.1038/s41598-017-05425-7}
#'
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' phases <- hilbert_phases(ts)
#' dim(phases)  # 10 x 200
hilbert_phases <- function(timeseries) {
  N    <- nrow(timeseries)
  Tmax <- ncol(timeseries)
  phases <- matrix(0, nrow = N, ncol = Tmax)
  for (i in seq_len(N)) {
    x_c <- timeseries[i, ] - mean(timeseries[i, ])
    xan <- gsignal::hilbert(x_c)
    phases[i, ] <- Arg(xan)
  }
  phases
}
