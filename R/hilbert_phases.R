#' Hilbert transform phase extraction
#'
#' Compute the instantaneous phase time series for all parcels/voxels via
#' the analytic signal (Hilbert transform). Each parcel timeseries is
#' demeaned before transformation.
#'
#' The analytic signal is computed by `.hilbert_r()`, a bespoke port of
#' `scipy.signal.hilbert()` that uses only base-R `fft()` and requires no
#' external packages.
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
  N    <- nrow(timeseries)
  Tmax <- ncol(timeseries)
  phases <- matrix(0, nrow = N, ncol = Tmax)
  for (i in seq_len(N)) {
    x_c        <- timeseries[i, ] - mean(timeseries[i, ])
    phases[i, ] <- Arg(.hilbert_r(x_c))
  }
  phases
}

# ── Internal helpers (not exported) ──────────────────────────────────────────

# Analytic signal via FFT. Equivalent to scipy.signal.hilbert() and
# gsignal::hilbert(). Returns a complex vector of the same length as x.
#
# Algorithm: FFT -> zero the negative-frequency components (multiply by h)
# -> IFFT. The one-sided spectrum multiplier h is:
#   h[1]           = 1  (DC)
#   h[2 : N/2]     = 2  (positive freqs, doubled to preserve energy)
#   h[N/2 + 1]     = 1  (Nyquist, N even only)
#   h[N/2 + 2 : N] = 0  (negative freqs, zeroed)
# For odd N the Nyquist bin does not exist; positive freqs run to (N+1)/2.
.hilbert_r <- function(x) {
  n  <- length(x)
  Xf <- fft(x)
  h  <- numeric(n)
  if (n %% 2L == 0L) {
    h[1L]           <- 1.0           # DC
    h[n %/% 2L + 1L] <- 1.0          # Nyquist
    h[2L:(n %/% 2L)] <- 2.0          # positive frequencies
  } else {
    h[1L]                    <- 1.0   # DC
    h[2L:((n + 1L) %/% 2L)] <- 2.0   # positive frequencies
  }
  fft(Xf * h, inverse = TRUE) / n
}
