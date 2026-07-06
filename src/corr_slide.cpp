#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

//' Sliding-window Pearson correlation (C++ backend)
//'
//' Internal workhorse called by [corr_slide()]. Computes Pearson correlation
//' matrices directly on the input matrix without allocating per-window
//' submatrices or calling R's `cor()`.
//'
//' Loop ordering is chosen for cache efficiency on the column-major R matrix:
//' the innermost loop over `t` (timepoints) within a window iterates over
//' consecutive memory locations in each column of `timeseries`, so the
//' cross-product accumulation step is sequential in memory.
//'
//' @param timeseries NumericMatrix \[N x Tmax\].
//' @param window Integer window size in timepoints.
//' @param step Integer step between window onsets.
//'
//' @return List with:
//'   \item{corr_mats}{NumericVector with dim \[N, N, n_windows\].}
//'   \item{idx}{IntegerVector of 1-indexed window onset positions.}
//'
//' @keywords internal
// [[Rcpp::export]]
List corr_slide_cpp(NumericMatrix timeseries, int window, int step) {
  const int N    = timeseries.nrow();
  const int Tmax = timeseries.ncol();

  // ── Window onset indices (1-indexed, matching R's seq()) ─────────────────
  std::vector<int> idx;
  for (int i = 1; i + window - 1 <= Tmax; i += step)
    idx.push_back(i);
  const int n_windows = static_cast<int>(idx.size());

  // ── Output array N x N x n_windows (column-major) ────────────────────────
  NumericVector corr_mats(N * N * n_windows);

  // ── Per-window working storage (reused across windows) ───────────────────
  std::vector<double> means(N);
  std::vector<double> sds(N);
  std::vector<double> cp(N * N, 0.0);   // cross-product accumulator (dense)

  for (int w = 0; w < n_windows; w++) {
    const int t0    = idx[w] - 1;        // 0-indexed window start
    const int t_off = w * N * N;

    // ── Step 1: means ────────────────────────────────────────────────────────
    for (int i = 0; i < N; i++) {
      double s = 0.0;
      for (int t = t0; t < t0 + window; t++) s += timeseries(i, t);
      means[i] = s / window;
    }

    // ── Step 2: cross-products (t outer = column-major sequential access) ────
    // cp[i*N + j] accumulates sum_t (x[i,t]-mi)(x[j,t]-mj), upper triangle.
    // Zero first (reuse allocation).
    std::fill(cp.begin(), cp.end(), 0.0);
    for (int t = t0; t < t0 + window; t++) {
      for (int i = 0; i < N; i++) {
        const double xi = timeseries(i, t) - means[i];
        for (int j = i; j < N; j++) {
          cp[i * N + j] += xi * (timeseries(j, t) - means[j]);
        }
      }
    }

    // ── Step 3: standard deviations (diagonal of cp, scaled) ─────────────────
    for (int i = 0; i < N; i++)
      sds[i] = std::sqrt(cp[i * N + i] / (window - 1));

    // ── Step 4: correlation matrix ────────────────────────────────────────────
    for (int i = 0; i < N; i++) {
      corr_mats[t_off + i + i * N] = 1.0;
      for (int j = i + 1; j < N; j++) {
        double c = (cp[i * N + j] / (window - 1)) / (sds[i] * sds[j]);
        corr_mats[t_off + i + j * N] = c;
        corr_mats[t_off + j + i * N] = c;
      }
    }
  }

  corr_mats.attr("dim") = IntegerVector::create(N, N, n_windows);

  IntegerVector r_idx(n_windows);
  for (int w = 0; w < n_windows; w++) r_idx[w] = idx[w];

  return List::create(Named("corr_mats") = corr_mats, Named("idx") = r_idx);
}
