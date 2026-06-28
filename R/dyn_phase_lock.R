#' Dynamic phase-locking matrix (dPL)
#'
#' Compute instantaneous phase-locking matrices from parcel-level phase time
#' series and extract leading eigenvectors via [get_leida()]. The first and
#' last 10 timepoints are discarded to avoid edge effects from the Hilbert
#' transform.
#'
#' @param phases Numeric matrix \[N × Tmax\]. Instantaneous phases in radians,
#'   as returned by [hilbert_phases()].
#'
#' @return A list with:
#'   \item{sync_conn}{Array \[N, N, Tmax-20\]. Instantaneous phase-locking
#'     matrices, one per (trimmed) timepoint.}
#'   \item{leida}{Matrix \[Tmax-20, N\]. Leading eigenvectors from [get_leida()].}
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
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10, ncol = 200)
#' phases <- hilbert_phases(ts)
#' res <- dyn_phase_lock(phases)
#' dim(res$sync_conn)  # 10 x 10 x 180
#' dim(res$leida)      # 180 x 10
dyn_phase_lock <- function(phases) {
  N     <- nrow(phases)
  Tmax  <- ncol(phases)
  T_idx <- seq(11L, Tmax - 10L)   # trim 10 timepoints from each end
  n_T   <- length(T_idx)
  sync_conn <- array(0, dim = c(N, N, n_T))
  for (t in seq_len(n_T)) {
    ph_t <- phases[, T_idx[t]]
    sync_conn[, , t] <- outer(ph_t, ph_t, function(a, b) cos(a - b))
  }
  leida <- get_leida(sync_conn)
  list(sync_conn = sync_conn, leida = leida)
}
