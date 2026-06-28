#' BOLD fMRI timeseries (200 parcels, 600 timepoints)
#'
#' A resting-state BOLD fMRI timeseries from the
#' [edge-ts](https://github.com/brain-networks/edge-ts) repository, parcellated
#' into 200 regions of interest. Used in dynR vignettes and examples to
#' demonstrate dynamic FC methods on realistic data.
#'
#' @format A numeric matrix with 200 rows (parcels) and 600 columns (timepoints).
#'
#' @source \url{https://github.com/brain-networks/edge-ts}
"ts"

#' Functional connectivity matrix (200 parcels)
#'
#' Static functional connectivity matrix derived from [ts] by computing the
#' full-timeseries Pearson correlation across all 200 parcel pairs. Serves as
#' a ground-truth reference for validating sliding-window and phase-based dynFC
#' methods when the window spans the full timeseries.
#'
#' @format A numeric matrix with 200 rows and 200 columns. Diagonal is 1;
#'   off-diagonal values are Pearson correlations in \[-1, 1\].
#'
#' @source Derived from [ts].
"fc"
