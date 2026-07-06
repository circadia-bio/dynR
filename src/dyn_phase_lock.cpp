#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Instantaneous phase-locking matrices (C++ backend)
//'
//' Internal workhorse called by [dyn_phase_lock()]. Iterates over the trimmed
//' timepoints (dropping the first and last 10) and fills an \[N, N, n_T\]
//' array with \eqn{\cos(\phi_i - \phi_j)} for every parcel pair.
//'
//' Exploits the symmetry \eqn{\cos(a-b) = \cos(b-a)}: only the upper triangle
//' is computed via \code{std::cos}; the diagonal is set to 1 directly, and
//' values are mirrored to the lower triangle. This halves the number of
//' trigonometric evaluations relative to a naive double loop.
//'
//' @param phases NumericMatrix \[N x Tmax\]. Instantaneous phases in radians,
//'   as returned by [hilbert_phases()].
//'
//' @return A NumericVector with a \code{dim} attribute \[N, N, n_T\] (i.e. a
//'   3-D array), where \code{n_T = Tmax - 20}.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector dyn_phase_lock_cpp(NumericMatrix phases) {
  const int N      = phases.nrow();
  const int Tmax   = phases.ncol();
  const int t_start = 10;           // 0-indexed: corresponds to R timepoint 11
  const int t_end   = Tmax - 11;    // 0-indexed: corresponds to R timepoint Tmax-10
  const int n_T     = t_end - t_start + 1;

  NumericVector sync_conn(N * N * n_T);

  for (int t = 0; t < n_T; t++) {
    const int t_idx = t_start + t;
    const int t_off = t * N * N;

    for (int i = 0; i < N; i++) {
      // Diagonal: cos(0) = 1
      sync_conn[t_off + i + i * N] = 1.0;

      // Upper triangle + mirror to lower
      for (int j = i + 1; j < N; j++) {
        double val = std::cos(phases(i, t_idx) - phases(j, t_idx));
        sync_conn[t_off + i + j * N] = val;
        sync_conn[t_off + j + i * N] = val;
      }
    }
  }

  sync_conn.attr("dim") = IntegerVector::create(N, N, n_T);
  return sync_conn;
}
