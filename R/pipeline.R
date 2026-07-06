# dynR colour palette — hard-coded so circadia is not a package dependency.
# Matches the pkgdown theme and both vignettes exactly.
# State colours matching both vignettes (up to 5 states; recycled if K > 5):
.state_cols <- function(k) {
  rep_len(c("#8D9FD7", "#9E3C30", "#754D71", "#341C5D", "#E19A8F"), k)
}

# Validate filter parameters and error with modality guidance if any are missing.
.check_filter_params <- function(flp, fhi, delt) {
  missing_args <- c(
    if (is.null(flp))  "flp",
    if (is.null(fhi))  "fhi",
    if (is.null(delt)) "delt"
  )
  if (length(missing_args) == 0L) return(invisible(NULL))

  stop(
    "filter = TRUE but the following arguments are not set: ",
    paste(missing_args, collapse = ", "), ".\n",
    "Set them explicitly for your recording modality, for example:\n",
    "  fMRI BOLD  : flp = 0.01, fhi = 0.1,  delt = <TR in seconds>\n",
    "  EEG alpha  : flp = 8,    fhi = 12,   delt = 1 / <sampling_rate>\n",
    "  EEG theta  : flp = 4,    fhi = 8,    delt = 1 / <sampling_rate>\n",
    "  LFP gamma  : flp = 30,   fhi = 80,   delt = 1 / <sampling_rate>\n",
    "Or pass filter = FALSE if the timeseries is already filtered.",
    call. = FALSE
  )
}

#' LEiDA pipeline
#'
#' Convenience wrapper running the full phase-based dynamic FC pipeline:
#' optional Butterworth bandpass filter -> Hilbert transform ->
#' dynamic phase-locking matrices + LEiDA eigenvectors -> Kuramoto order
#' parameter. Returns a structured `dynR_leida` object with `print()` and
#' `plot()` methods.
#'
#' dynR is modality-agnostic. `flp`, `fhi`, and `delt` have no defaults and
#' **must** be supplied when `filter = TRUE`. Set them for your recording
#' modality (see examples). Pass `filter = FALSE` if the timeseries is already
#' band-limited.
#'
#' @param timeseries Numeric matrix \[N x Tmax\]. Rows = channels/parcels,
#'   columns = timepoints.
#' @param flp Numeric. Low-pass cutoff frequency (Hz). No default — must be
#'   set explicitly when `filter = TRUE`.
#' @param fhi Numeric. High-pass cutoff frequency (Hz). No default — must be
#'   set explicitly when `filter = TRUE`.
#' @param delt Numeric. Sampling interval in seconds (`1 / sampling_rate`).
#'   No default — must be set explicitly when `filter = TRUE`.
#' @param order Integer. Butterworth filter order. Default `2L`.
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
#'   \item{N}{Integer. Number of channels/parcels.}
#'   \item{Tmax}{Integer. Number of timepoints.}
#'
#' @seealso `sw_pipeline()`, `dyn_phase_lock()`, `kuramoto()`, `plot.dynR_leida()`
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10)
#'
#' # filter = FALSE: timeseries already band-limited
#' res <- leida_pipeline(ts, filter = FALSE)
#' res
#'
#' # fMRI BOLD (TR = 2 s)
#' \dontrun{
#' res <- leida_pipeline(ts, flp = 0.01, fhi = 0.1, delt = 2)
#' }
#'
#' # EEG alpha band (250 Hz)
#' \dontrun{
#' res <- leida_pipeline(ts, flp = 8, fhi = 12, delt = 1 / 250)
#' }
leida_pipeline <- function(timeseries, flp = NULL, fhi = NULL, delt = NULL,
                            order = 2L, filter = TRUE) {
  if (filter) {
    .check_filter_params(flp, fhi, delt)
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
#' optional Butterworth bandpass filter -> sliding-window Pearson correlation
#' matrices -> edge-centric cofluctuations. Returns a structured `dynR_sw`
#' object with `print()` and `plot()` methods.
#'
#' dynR is modality-agnostic. `flp`, `fhi`, and `delt` have no defaults and
#' **must** be supplied when `filter = TRUE`. Pass `filter = FALSE` if the
#' timeseries is already band-limited.
#'
#' @param timeseries Numeric matrix \[N x Tmax\].
#' @param window Integer. Window size in timepoints.
#' @param step Integer. Step between window onsets. Default: `window`
#'   (non-overlapping windows).
#' @param flp Numeric. Low-pass cutoff (Hz). No default — required when
#'   `filter = TRUE`.
#' @param fhi Numeric. High-pass cutoff (Hz). No default — required when
#'   `filter = TRUE`.
#' @param delt Numeric. Sampling interval in seconds. No default — required
#'   when `filter = TRUE`.
#' @param order Integer. Butterworth filter order. Default `2L`.
#' @param filter Logical. Apply bandpass filter? Default `TRUE`.
#'
#' @return An object of class `dynR_sw` (a named list) with elements:
#'   \item{corr_mats}{Array \[N, N, n_windows\]. Sliding FC matrices.}
#'   \item{idx}{Integer vector. 1-indexed window onset positions.}
#'   \item{edge_ts}{Matrix \[n_edges, Tmax\]. Edge time series.}
#'   \item{rss}{Numeric vector \[Tmax\]. Root-sum-square cofluctuation.}
#'   \item{window}{Integer. Window size used.}
#'   \item{step}{Integer. Step size used.}
#'   \item{N}{Integer. Number of channels/parcels.}
#'   \item{Tmax}{Integer. Number of timepoints.}
#'
#' @seealso `leida_pipeline()`, `corr_slide()`, `cofluct()`, `plot.dynR_sw()`
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10)
#'
#' # filter = FALSE: timeseries already band-limited
#' res <- sw_pipeline(ts, window = 20, filter = FALSE)
#' res
#'
#' # fMRI BOLD (TR = 2 s)
#' \dontrun{
#' res <- sw_pipeline(ts, window = 20, flp = 0.01, fhi = 0.1, delt = 2)
#' }
sw_pipeline <- function(timeseries, window, step = NULL, flp = NULL, fhi = NULL,
                         delt = NULL, order = 2L, filter = TRUE) {
  if (is.null(step)) step <- window
  if (filter) {
    .check_filter_params(flp, fhi, delt)
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
