#' Leading eigenvector decomposition (LEiDA)
#'
#' Extract the leading eigenvector from a series of instantaneous
#' phase-locking matrices. The sign of each eigenvector is normalised so that
#' its sum is non-positive, following the LEiDA convention.
#'
#' The computation is accelerated by a compiled C++ / LAPACK backend
#' ([get_leida_cpp()]): LAPACK `dsyev` is called once per timepoint with a
#' single shared workspace allocation, and only the leading eigenvector
#' (last column of the ascending-order output) is retained.
#'
#' @param sync_conn Numeric array \[N, N, Tmax\]. Instantaneous phase-locking
#'   matrices, as returned by [dyn_phase_lock()].
#'
#' @return Numeric matrix \[Tmax x N\]. Leading eigenvectors, one per
#'   timepoint.
#'
#' @references
#' Cabral, J. et al. (2017). Cognitive performance in healthy older adults
#' relates to spontaneous switching between states of functional connectivity
#' during rest. *Scientific Reports*, 7(1), 5135.
#' \doi{10.1038/s41598-017-05425-7}
#'
#' Lord, L.-D. et al. (2019). Dynamical exploration of the repertoire of brain
#' networks at rest is modulated by psilocybin. *NeuroImage*, 199, 127-142.
#' \doi{10.1016/j.neuroimage.2019.05.060}
#'
#' @export
get_leida <- function(sync_conn) {
  get_leida_cpp(sync_conn)
}
