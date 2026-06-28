#' Leading eigenvector decomposition (LEiDA)
#'
#' Extract the leading eigenvector from a series of instantaneous
#' phase-locking matrices. The sign of each eigenvector is normalised so that
#' its sum is non-positive, following the LEiDA convention.
#'
#' @param sync_conn Numeric array \[N, N, Tmax\]. Instantaneous phase-locking
#'   matrices, as returned by [dyn_phase_lock()].
#'
#' @return Numeric matrix \[Tmax × N\]. Leading eigenvectors, one per
#'   timepoint.
#'
#' @references
#' Cabral, J. et al. (2017). Cognitive performance in healthy older adults
#' relates to spontaneous switching between states of functional connectivity
#' during rest. *Scientific Reports*, 7(1), 5135.
#' \doi{10.1038/s41598-017-05425-7}
#'
#' Lord, L.-D. et al. (2019). Dynamical exploration of the repertoire of brain
#' networks at rest is modulated by psilocybin. *NeuroImage*, 199, 127–142.
#' \doi{10.1016/j.neuroimage.2019.05.060}
#'
#' @export
get_leida <- function(sync_conn) {
  t_points <- dim(sync_conn)[3]
  N        <- dim(sync_conn)[1]
  leida    <- matrix(0, nrow = t_points, ncol = N)
  for (i in seq_len(t_points)) {
    ev <- eigen(sync_conn[, , i], symmetric = TRUE)
    v1 <- ev$vectors[, 1]   # eigen() orders decreasing — [,1] is the leading eigenvector
    if (sum(v1) > 0) v1 <- -v1
    leida[i, ] <- v1
  }
  leida
}
