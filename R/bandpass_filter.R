#' Butterworth bandpass filter
#'
#' Design and apply a zero-phase Butterworth bandpass filter to a signal.
#' Wraps [gsignal::butter()] and [gsignal::filtfilt()] for convenience.
#'
#' @param x Numeric vector. Signal to be filtered.
#' @param flp Numeric. Low-pass cutoff frequency (Hz).
#' @param fhi Numeric. High-pass cutoff frequency (Hz).
#' @param delt Numeric. Sampling interval in seconds (i.e. TR for fMRI).
#' @param order Integer. Filter order. Default is 2.
#'
#' @return Numeric vector. Zero-phase filtered signal, same length as `x`.
#'
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(200)
#' x_filt <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 2, order = 2)
bandpass_filter <- function(x, flp, fhi, delt, order = 2) {
  fnq  <- 1 / (2 * delt)
  low  <- flp / fnq
  high <- fhi / fnq
  filt <- gsignal::butter(order, c(low, high), type = "pass")
  as.numeric(gsignal::filtfilt(filt, x))
}
