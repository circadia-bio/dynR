#' Dynamic phase-locking matrix (dPL)
#'
#' Compute instantaneous phase-locking matrices from parcel-level phase time
#' series and extract leading eigenvectors via [get_leida()]. The first and
#' last 10 timepoints are discarded to avoid edge effects from the Hilbert
#' transform.
#'
#' The phase-locking computation is accelerated by a compiled C++ backend
#' (`dyn_phase_lock_cpp()`): only the upper triangle is evaluated with
#' `cos(phi_i - phi_j)` and mirrored to the lower, halving trigonometric
#' operations relative to a naive double loop.
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
#' networks at rest is modulated by psilocybin. *NeuroImage*, 199, 127-142.
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
  sync_conn <- dyn_phase_lock_cpp(phases)
  leida     <- get_leida(sync_conn)
  list(sync_conn = sync_conn, leida = leida)
}
