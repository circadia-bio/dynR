# dynR colour palette — hard-coded so circadia is not a package dependency.
# Matches the pkgdown theme and both vignettes exactly.
# State colours matching both vignettes (up to 5 states; recycled if K > 5):
.state_cols <- function(k) {
  rep_len(c("#8D9FD7", "#9E3C30", "#754D71", "#341C5D", "#E19A8F"), k)
}

#' LEiDA pipeline
#'
#' Convenience wrapper running the full phase-based dynamic FC pipeline:
#' optional Butterworth bandpass filter → Hilbert transform →
#' dynamic phase-locking matrices + LEiDA eigenvectors → Kuramoto order
#' parameter. Returns a structured `dynR_leida` object with `print()` and
#' `plot()` methods.
#'
#' @param timeseries Numeric matrix \[N x Tmax\]. BOLD signal with N parcels
#'   as rows and Tmax timepoints as columns.
#' @param flp Numeric. Low-pass cutoff frequency (Hz). Default 0.01.
#' @param fhi Numeric. High-pass cutoff frequency (Hz). Default 0.1.
#' @param delt Numeric. Sampling interval in seconds (TR). Default 2.
#' @param order Integer. Butterworth filter order. Default 2.
#' @param filter Logical. Apply bandpass filter before phase extraction?
#'   Default `TRUE`. Set to `FALSE` if the timeseries is already filtered.
#'
#' @return An object of class `dynR_leida` (a named list) with elements:
#'   \item{leida}{Matrix \[Tmax-20, N\]. LEiDA eigenvectors.}
#'   \item{sync_conn}{Array \[N, N, Tmax-20\]. Instantaneous phase-locking matrices.}
#'   \item{phases}{Matrix \[N, Tmax\]. Instantaneous phases in radians.}
#'   \item{synchrony}{Numeric vector \[Tmax-20\]. Kuramoto R(t).}
#'   \item{metastability}{Numeric. Standard deviation of synchrony.}
#'   \item{entropy}{Numeric. Shannon entropy of synchrony series (bits).}
#'   \item{N}{Integer. Number of parcels.}
#'   \item{Tmax}{Integer. Number of timepoints.}
#'
#' @seealso `sw_pipeline()`, `dyn_phase_lock()`, `kuramoto()`, `plot.dynR_leida()`
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10)
#' res <- leida_pipeline(ts, filter = FALSE)
#' res
leida_pipeline <- function(timeseries, flp = 0.01, fhi = 0.1, delt = 2,
                            order = 2L, filter = TRUE) {
  if (filter) {
    timeseries <- t(apply(timeseries, 1L, bandpass_filter,
                          flp = flp, fhi = fhi, delt = delt, order = order))
  }
  phases <- hilbert_phases(timeseries)
  dpl    <- dyn_phase_lock(phases)
  kop    <- kuramoto(phases)

  structure(
    list(
      leida         = dpl$leida,
      sync_conn     = dpl$sync_conn,
      phases        = phases,
      synchrony     = kop$synchrony,
      metastability = kop$metastability,
      entropy       = kop$entropy,
      N             = nrow(timeseries),
      Tmax          = ncol(timeseries)
    ),
    class = "dynR_leida"
  )
}

#' Sliding-window pipeline
#'
#' Convenience wrapper running the full correlation-based dynamic FC pipeline:
#' optional Butterworth bandpass filter → sliding-window Pearson correlation
#' matrices → edge-centric cofluctuations. Returns a structured `dynR_sw`
#' object with `print()` and `plot()` methods.
#'
#' @param timeseries Numeric matrix \[N x Tmax\].
#' @param window Integer. Window size in timepoints.
#' @param step Integer. Step between window onsets. Default: `window`
#'   (non-overlapping windows).
#' @param flp,fhi,delt,order Filter parameters passed to [bandpass_filter()].
#' @param filter Logical. Apply bandpass filter? Default `TRUE`.
#'
#' @return An object of class `dynR_sw` (a named list) with elements:
#'   \item{corr_mats}{Array \[N, N, n_windows\]. Sliding FC matrices.}
#'   \item{idx}{Integer vector. 1-indexed window onset positions.}
#'   \item{edge_ts}{Matrix \[n_edges, Tmax\]. Edge time series.}
#'   \item{rss}{Numeric vector \[Tmax\]. Root-sum-square cofluctuation.}
#'   \item{window}{Integer. Window size used.}
#'   \item{step}{Integer. Step size used.}
#'   \item{N}{Integer. Number of parcels.}
#'   \item{Tmax}{Integer. Number of timepoints.}
#'
#' @seealso `leida_pipeline()`, `corr_slide()`, `cofluct()`, `plot.dynR_sw()`
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10)
#' res <- sw_pipeline(ts, window = 20, filter = FALSE)
#' res
sw_pipeline <- function(timeseries, window, step = NULL, flp = 0.01, fhi = 0.1,
                         delt = 2, order = 2L, filter = TRUE) {
  if (is.null(step)) step <- window
  if (filter) {
    timeseries <- t(apply(timeseries, 1L, bandpass_filter,
                          flp = flp, fhi = fhi, delt = delt, order = order))
  }
  sw <- corr_slide(timeseries, window = window, step = step)
  ec <- cofluct(timeseries)

  structure(
    list(
      corr_mats = sw$corr_mats,
      idx       = sw$idx,
      edge_ts   = ec$edge_ts,
      rss       = ec$rss,
      window    = window,
      step      = step,
      N         = nrow(timeseries),
      Tmax      = ncol(timeseries)
    ),
    class = "dynR_sw"
  )
}

#' @export
print.dynR_leida <- function(x, ...) {
  cat("<dynR_leida>\n")
  cat("  Parcels:      ", x$N, "\n")
  cat("  Timepoints:   ", x$Tmax, "\n")
  cat("  LEiDA frames: ", nrow(x$leida), "\n")
  cat("  Metastability:", round(x$metastability, 4), "\n")
  cat("  Entropy:      ", round(x$entropy, 4), "bits\n")
  invisible(x)
}

#' @export
print.dynR_sw <- function(x, ...) {
  cat("<dynR_sw>\n")
  cat("  Parcels:    ", x$N, "\n")
  cat("  Timepoints: ", x$Tmax, "\n")
  cat("  Window:     ", x$window, "timepoints\n")
  cat("  Step:       ", x$step, "timepoints\n")
  cat("  N windows:  ", length(x$idx), "\n")
  cat("  N edges:    ", nrow(x$edge_ts), "\n")
  invisible(x)
}
