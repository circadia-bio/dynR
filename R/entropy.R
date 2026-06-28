#' Shannon entropy
#'
#' Estimate Shannon entropy from a numeric vector, with optional
#' bit-depth discretisation (as used internally by [kuramoto()]).
#'
#' @param x Numeric vector.
#' @param base Numeric. Logarithm base. Default is 2 (entropy in bits).
#' @param n_bits Integer or `NULL`. If supplied, values are scaled by
#'   \eqn{2^{n\_bits}} and rounded before computing entropy. Default is `NULL`.
#'
#' @return Numeric scalar. Shannon entropy.
#'
#' @export
#' @examples
#' x <- sample(1:4, 100, replace = TRUE)
#' shannon_entropy(x)
shannon_entropy <- function(x, base = 2, n_bits = NULL) {
  if (!is.null(n_bits)) x <- round(x * 2^n_bits)
  probs <- tabulate(factor(x)) / length(x)
  probs <- probs[probs > 0]
  -sum(probs * log(probs, base = base))
}
