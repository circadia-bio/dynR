#' Hilbert transform phase extraction
#'
#' Compute the instantaneous phase time series for all parcels/voxels via
#' the analytic signal (Hilbert transform). Each parcel timeseries is
#' demeaned before transformation.
#'
#' All N parcels are processed in two vectorised calls to [mvfft()] rather
#' than an N-iteration R loop over [fft()]: the timeseries matrix is transposed
#' to a Tmax x N layout so `mvfft` applies the FFT to each parcel column
#' simultaneously, the Hilbert multiplier is broadcast across columns, and
#' the inverse FFT recovers the analytic signal for all parcels at once.
#'
#' `.hilbert_r()` is retained as an internal single-vector reference
#' (scipy-parity validated).
#'
#' @param timeseries Numeric matrix \[N x Tmax\]. BOLD signal with N parcels
#'   as rows and Tmax timepoints as columns.
#'
#' @return Numeric matrix \[N x Tmax\]. Instantaneous phases in radians.
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
  Tmax <- ncol(timeseries)

  # Demean each parcel (row) in one vectorised step
  ts_zero <- timeseries - rowMeans(timeseries)

  # mvfft expects columns -> transpose to Tmax x N
  Xf <- mvfft(t(ts_zero))

  # Hilbert one-sided multiplier (same for every parcel)
  h <- numeric(Tmax)
  h[1L] <- 1.0
  if (Tmax %% 2L == 0L) {
    h[Tmax %/% 2L + 1L] <- 1.0
    h[2L:(Tmax %/% 2L)] <- 2.0
  } else {
    h[2L:((Tmax + 1L) %/% 2L)] <- 2.0
  }

  # h is recycled column-by-column over the Tmax x N matrix: each column
  # (parcel) gets multiplied by the same h vector
  Xan <- mvfft(Xf * h, inverse = TRUE) / Tmax

  # Phases: Arg of analytic signal, transposed back to N x Tmax
  t(Arg(Xan))
}

# ── Internal reference implementation (not exported) ─────────────────────────

# Single-vector analytic signal via FFT. Equivalent to scipy.signal.hilbert().
# Retained for parity validation; hilbert_phases() uses the vectorised path.
.hilbert_r <- function(x) {
  n  <- length(x)
  Xf <- fft(x)
  h  <- numeric(n)
  if (n %% 2L == 0L) {
    h[1L]             <- 1.0
    h[n %/% 2L + 1L]  <- 1.0
    h[2L:(n %/% 2L)]  <- 2.0
  } else {
    h[1L]                    <- 1.0
    h[2L:((n + 1L) %/% 2L)] <- 2.0
  }
  fft(Xf * h, inverse = TRUE) / n
}
