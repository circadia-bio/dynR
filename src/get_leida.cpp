#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

//' Leading eigenvector decomposition (LEiDA) -- C++ / LAPACK backend
//'
//' Internal workhorse called by [get_leida()]. For each phase-locking matrix
//' slice in `sync_conn`, calls LAPACK `dsyev` (symmetric eigendecomposition)
//' and returns only the leading eigenvector.
//'
//' Differences from the R implementation:
//' - `dsyev` returns eigenvalues in **ascending** order, so the leading
//'   eigenvector is the **last** column of the output (opposite of R's
//'   `eigen()` which is descending).
//' - The LAPACK workspace is allocated once and reused across all timepoints.
//' - Sign convention preserved: row sum forced to be non-positive.
//'
//' @param sync_conn NumericVector with `dim` attribute \[N, N, t_points\].
//'
//' @return NumericMatrix \[t_points x N\]. Leading eigenvectors, one per row.
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix get_leida_cpp(NumericVector sync_conn) {
  IntegerVector d    = sync_conn.attr("dim");
  const int N        = d[0];
  const int t_points = d[2];

  NumericMatrix leida(t_points, N);

  // dsyev arguments
  char jobz = 'V';   // compute eigenvalues AND eigenvectors
  char uplo = 'U';   // reference the upper triangle
  int  n    = N;
  int  lda  = N;
  int  info = 0;

  std::vector<double> evals(N);
  std::vector<double> A(N * N);

  // ── Workspace query ───────────────────────────────────────────────────────
  // Call dsyev with LWORK = -1; it returns the optimal workspace size in
  // WORK[0] without doing any real computation.
  {
    int    lwork_q = -1;
    double work_q  = 0.0;
    A.assign(N * N, 0.0);
    F77_CALL(dsyev)(&jobz, &uplo, &n, A.data(), &lda,
                    evals.data(), &work_q, &lwork_q, &info FCONE FCONE);
    int lwork = static_cast<int>(work_q);
    if (lwork < 1) lwork = 3 * N + 64;  // safe fallback

    std::vector<double> work(lwork);

    // ── Main loop: one dsyev call per timepoint ───────────────────────────
    for (int t = 0; t < t_points; t++) {
      // Copy slice t (column-major N x N) from sync_conn.
      // R array [N, N, T] stores element [i,j,t] at i + j*N + t*N*N.
      const double* src = &sync_conn[t * N * N];
      std::copy(src, src + N * N, A.begin());

      // dsyev overwrites A with eigenvectors; evals filled ascending.
      info = 0;
      F77_CALL(dsyev)(&jobz, &uplo, &n, A.data(), &lda,
                      evals.data(), work.data(), &lwork, &info FCONE FCONE);

      if (info != 0) {
        Rcpp::warning("dsyev failed at timepoint %d (info = %d)", t + 1, info);
      }

      // Leading eigenvector = last column of A (ascending eigenvalue order)
      const double* v1 = &A[(N - 1) * N];

      // Sign convention: force row sum <= 0 (LEiDA standard)
      double s = 0.0;
      for (int i = 0; i < N; i++) s += v1[i];
      const double sgn = (s > 0.0) ? -1.0 : 1.0;

      for (int i = 0; i < N; i++) leida(t, i) = sgn * v1[i];
    }
  }

  return leida;
}
