#' Euclidean distance between consecutive points
#'
#' Compute the Euclidean distance between each consecutive pair of rows in a
#' matrix — typically a trajectory through PC space.
#'
#' @param x Numeric matrix \[n_points × n_dims\].
#'
#' @return Numeric vector of length `nrow(x)`. The first element is always 0
#'   (no previous point).
#'
#' @export
#' @examples
#' set.seed(1)
#' pcs <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
#' d <- do_euclid(pcs)
#' length(d)  # 50
do_euclid <- function(x) {
  n <- nrow(x)
  d <- numeric(n)
  for (i in seq(2L, n)) {
    d[i] <- sqrt(sum((x[i, ] - x[i - 1L, ])^2))
  }
  d
}
