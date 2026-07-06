#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Kuramoto order parameter time series (C++ backend)
//'
//' Internal workhorse called by [kuramoto()]. For each trimmed timepoint,
//' computes the global Kuramoto order parameter R(t):
//'
//' \deqn{R(t) = \frac{1}{N} \left| \sum_{j=1}^{N} e^{i\phi_j(t)} \right|}
//'
//' The complex exponential is decomposed as \eqn{\cos\phi + i\sin\phi};
//' the modulus is \eqn{\sqrt{(\sum\cos)^2 + (\sum\sin)^2} / N}, avoiding
//' any R complex-number allocation.
//'
//' The first and last 10 timepoints are discarded (matching [dyn_phase_lock()]).
//'
//' @param phases NumericMatrix \[N x Tmax\]. Instantaneous phases in radians.
//'
//' @return NumericVector \[Tmax-20\]. Kuramoto order parameter per timepoint.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector kuramoto_sync_cpp(NumericMatrix phases) {
  const int N       = phases.nrow();
  const int Tmax    = phases.ncol();
  const int t_start = 10;           // 0-indexed: R's timepoint 11
  const int t_end   = Tmax - 11;    // 0-indexed: R's timepoint Tmax-10
  const int n_T     = t_end - t_start + 1;

  NumericVector sync(n_T);

  for (int t = 0; t < n_T; t++) {
    const int t_idx = t_start + t;
    double sc = 0.0;   // sum of cos(phi)
    double ss = 0.0;   // sum of sin(phi)
    for (int j = 0; j < N; j++) {
      sc += std::cos(phases(j, t_idx));
      ss += std::sin(phases(j, t_idx));
    }
    sync[t] = std::sqrt(sc * sc + ss * ss) / N;
  }

  return sync;
}
