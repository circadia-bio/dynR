#' Kuramoto order parameter and metastability
#'
#' Compute the global Kuramoto order parameter time series, the metastability
#' index (standard deviation of the order parameter), and Shannon entropy of
#' synchrony from parcel-level instantaneous phase time series.
#' The first and last 10 timepoints are discarded (matching [dyn_phase_lock()]).
#'
#' @param phases Numeric matrix \[N × Tmax\]. Instantaneous phases in radians,
#'   as returned by [hilbert_phases()].
#' @param base Numeric. Logarithm base for Shannon entropy. Default is 2.
#' @param n_bits Integer or `NULL`. Bit depth for discretising the synchrony
#'   series before entropy estimation. Default is 8.
#'
#' @return A list with:
#'   \item{metastability}{Numeric scalar. Standard deviation of the Kuramoto
#'     order parameter.}
#'   \item{synchrony}{Numeric vector \[Tmax-20\]. Kuramoto order parameter
#'     at each (trimmed) timepoint.}
#'   \item{entropy}{Numeric scalar. Shannon entropy of the synchrony series.}
#'
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' phases <- hilbert_phases(ts)
#' res <- kuramoto(phases)
#' res$metastability
kuramoto <- function(phases, base = 2, n_bits = 8L) {
  N     <- nrow(phases)
  Tmax  <- ncol(phases)
  T_idx <- seq(11L, Tmax - 10L)
  sync  <- vapply(T_idx, function(t) {
    ku <- sum(exp(1i * phases[, t])) / N
    Mod(ku)
  }, numeric(1))
  list(
    metastability = sd(sync),
    synchrony     = sync,
    entropy       = shannon_entropy(sync, base = base, n_bits = n_bits)
  )
}
