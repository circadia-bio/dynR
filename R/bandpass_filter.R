#' Butterworth bandpass filter
#'
#' Design and apply a zero-phase Butterworth bandpass filter to a signal,
#' using scipy-compatible steady-state initial conditions (equivalent to
#' `scipy.signal.sosfiltfilt` / `scipy.signal.filtfilt` with
#' `padtype = "odd"`).
#'
#' Two internal code paths are used depending on `order`:
#'
#' * **`order <= 3`** (typical fMRI use): direct-form b/a coefficients with
#'   the companion-matrix initial-condition method. scipy `filtfilt`-parity
#'   validated to `< 1e-9` for standard fMRI parameters.
#'
#' * **`order >= 4`**: second-order sections (SOS). The companion matrix for
#'   the full-order b/a representation becomes ill-conditioned for high-order
#'   filters with poles close to z = 1 (very low lower cutoff). SOS avoids
#'   this by decomposing the filter into biquad sections, each with a 2x2
#'   companion matrix that is always well-conditioned.
#'
#' Filter coefficients are computed by `.butter_bandpass_zpk()`, a bespoke
#' port of `scipy.signal.butter()` that requires no external packages.
#'
#' @param x Numeric vector. Signal to be filtered.
#' @param flp Numeric. Low-pass cutoff frequency (Hz).
#' @param fhi Numeric. High-pass cutoff frequency (Hz).
#' @param delt Numeric. Sampling interval in seconds (i.e. TR for fMRI).
#' @param order Integer. Filter order. Default is 2.
#'
#' @return Numeric vector. Zero-phase filtered signal, same length as `x`.
#'
#' @export
#' @examples
#' set.seed(1)
#' x <- rnorm(200)
#' x_filt <- bandpass_filter(x, flp = 0.01, fhi = 0.1, delt = 2, order = 2)
bandpass_filter <- function(x, flp, fhi, delt, order = 2) {
  fnq <- 1.0 / (2.0 * delt)
  low <- flp / fnq
  high <- fhi / fnq
  if (order <= 3L) {
    # b/a path: well-conditioned for low orders; scipy filtfilt-parity.
    filt <- .butter_bandpass(order, c(low, high))
    .filtfilt_zi(filt$b, filt$a, as.numeric(x))
  } else {
    # SOS path: numerically stable for high-order filters with poles near z=1.
    zpk <- .butter_bandpass_zpk(order, c(low, high))
    sos <- .zpk2sos_bp(zpk$z, zpk$p, zpk$k, order)
    .sosfiltfilt(sos, as.numeric(x))
  }
}

# ── Internal helpers (not exported) ──────────────────────────────────────────

# Zeros, poles, gain (ZPK) of a digital Butterworth bandpass filter.
# Shared by the b/a path (.butter_bandpass) and the SOS path (.zpk2sos_bp).
# Equivalent to scipy.signal.butter(nord, Wn, btype='bandpass', output='zpk').
#
# Pipeline: analog Butterworth LP prototype (butterap)
#           -> LP-to-BP transformation   (lp2bp_zpk)
#           -> bilinear transform         (bilinear_zpk, fs = 2)
#
# Wn: 2-element vector in [0, 1] where 1 = Nyquist.
.butter_bandpass_zpk <- function(nord, Wn) {
  fs  <- 2.0
  fs2 <- 4.0

  # Pre-warp (scipy: 2*fs*tan(pi*Wn/fs))
  warped <- fs2 * tan(pi * Wn / fs)
  bw <- warped[2L] - warped[1L]
  w0 <- sqrt(warped[1L] * warped[2L])

  # Analog Butterworth LP prototype (butterap)
  m_seq <- seq(-nord + 1L, nord - 1L, by = 2L)
  p_a   <- -exp(1i * pi * m_seq / (2.0 * nord))

  # LP -> BP transformation (lp2bp_zpk), no LP zeros
  hbw  <- as.complex(bw / 2.0) * p_a
  disc <- sqrt(hbw^2 - as.complex(w0^2))
  p_bp <- c(hbw + disc, hbw - disc)
  z_bp <- complex(real = numeric(nord))
  k_bp <- bw^nord

  # Bilinear transform (bilinear_zpk, fs = 2)
  # zeros at s=0 -> z=+1; zeros at infinity -> z=-1
  z_d <- c(rep(1 + 0i, nord), rep(-1 + 0i, nord))
  p_d <- (1.0 + p_bp / fs2) / (1.0 - p_bp / fs2)
  k_d <- k_bp * Re(prod(fs2 - z_bp) / prod(fs2 - p_bp))

  list(z = z_d, p = p_d, k = k_d)
}

# b/a coefficients from ZPK (used for order <= 3).
# Equivalent to scipy.signal.butter(order, Wn, btype='bandpass', output='ba').
.butter_bandpass <- function(order, Wn) {
  zpk <- .butter_bandpass_zpk(order, Wn)
  b   <- Re(.poly_from_roots(zpk$z)) * zpk$k
  a   <- Re(.poly_from_roots(zpk$p))
  a   <- a / a[1L]
  list(b = b, a = a)
}

# ── SOS path (order >= 4) ─────────────────────────────────────────────────────

# Convert Butterworth bandpass ZPK to SOS matrix (nord x 6).
# Butterworth bandpass structure:
#   zeros: nord at z = +1 and nord at z = -1 (all real, trivially paired as +-1)
#   poles: 2*nord complex values in nord conjugate pairs
# Each section: one (z=+1, z=-1) zero pair + one conjugate pole pair.
# Row format: [b0, b1, b2, 1, a1, a2].  Gain distributed onto the first row.
.zpk2sos_bp <- function(z, p, k, nord) {
  # Select one pole from each conjugate pair: sort all poles by Arg in [0, 2pi)
  # and take the first half (those with Arg in [0, pi)).
  angles <- Arg(p) %% (2 * pi)
  p_up   <- p[order(angles)][seq_len(nord)]

  sos <- matrix(0.0, nrow = nord, ncol = 6L)
  for (s in seq_len(nord)) {
    ps <- p_up[s]
    # Numerator: (x-1)(x+1) = x^2 - 1
    sos[s, 1L:3L] <- c(1.0, 0.0, -1.0)
    # Denominator: (x-ps)(x-conj(ps)) = x^2 - 2Re(ps)x + |ps|^2
    sos[s, 4L:6L] <- c(1.0, -2.0 * Re(ps), Re(ps)^2 + Im(ps)^2)
  }
  # Distribute gain onto the first section
  sos[1L, 1L:3L] <- sos[1L, 1L:3L] * k
  sos
}

# scipy.signal.sosfilt_zi equivalent.
# Computes steady-state initial conditions for all SOS sections with each
# section's zi scaled by the cumulative DC gain of all preceding sections.
# Returns an n_sections x 2 matrix (shape matches scipy's (n_sections, 2)).
#
# Port of scipy.signal.sosfilt_zi:
#   scale = 1.0
#   for section in range(n_sections):
#       zi[section] = scale * lfilter_zi(b, a)
#       scale *= sum(b) / sum(a)   # H(z=1): DC gain of this section
.sos_lfilter_zi <- function(sos) {
  ns    <- nrow(sos)
  zi    <- matrix(0.0, nrow = ns, ncol = 2L)
  scale <- 1.0
  for (s in seq_len(ns)) {
    b       <- sos[s, 1L:3L]
    a       <- sos[s, 4L:6L]
    zi[s, ] <- scale * .lfilter_zi(b, a)
    scale   <- scale * sum(b) / sum(a)
  }
  zi
}

# Cascade all SOS sections through x with per-section initial states zi.
# zi: n_sections x 2 matrix.  Port of scipy.signal.sosfilt (1-D).
# Output of each section becomes input to the next.
.sos_filt <- function(sos, x, zi) {
  ns <- nrow(sos)
  for (s in seq_len(ns)) {
    b       <- sos[s, 1L:3L]
    a       <- sos[s, 4L:6L]
    res     <- .lfilter(b, a, x, zi[s, ])
    x       <- res$y
    zi[s, ] <- res$zf
  }
  list(y = x, zf = zi)
}

# scipy.signal.sosfiltfilt equivalent with odd-reflection padding.
#
# padlen formula (direct port from scipy.signal.sosfiltfilt):
#   ntaps = 2*n_sections + 1 - min(sum(sos[:,2]==0), sum(sos[:,5]==0))
#   padlen = 3 * ntaps
# (The subtraction compensates for poles/zeros at the origin in odd-order filters.)
#
# zi (from .sos_lfilter_zi) is computed once then scaled by:
#   x[1]   for the forward pass
#   y[end] for the backward pass
# exactly mirroring scipy's sosfiltfilt.
.sosfiltfilt <- function(sos, x) {
  ns <- nrow(sos)

  # scipy padlen formula; columns 3 and 6 in R = scipy's 0-indexed columns 2 and 5
  n_b2 <- sum(sos[, 3L] == 0.0)
  n_a2 <- sum(sos[, 6L] == 0.0)
  pad  <- 3L * (2L * ns + 1L - min(n_b2, n_a2))

  ext <- .odd_ext(x, pad)

  # Initial conditions for all sections (DC-gain scaled)
  zi <- .sos_lfilter_zi(sos)

  # Forward pass: zi scaled by first sample of padded signal
  x0      <- ext[1L]
  res_fwd <- .sos_filt(sos, ext, zi * x0)

  # Backward pass: reverse signal, zi scaled by last sample of forward output
  y0      <- res_fwd$y[length(res_fwd$y)]
  res_bwd <- .sos_filt(sos, rev(res_fwd$y), zi * y0)

  # Reverse and strip padding
  y <- rev(res_bwd$y)
  y[seq(pad + 1L, pad + length(x))]
}

# ── b/a path helpers (order <= 3) ────────────────────────────────────────────

# Monic polynomial with the given roots, expanded via conjugate-pair real
# quadratics for numerical stability.
#
# Expanding a high-degree polynomial from complex linear factors one at a time
# introduces large cancellations that grow with filter order. For each pair of
# complex-conjugate roots (r, conj(r)), we instead multiply by the real
# quadratic  x^2 - 2*Re(r)*x + |r|^2  directly, keeping all arithmetic real
# and avoiding cancellation. Real roots are handled as linear factors.
#
# NOTE: R's convolve() computes cross-correlation by default (conj=TRUE), not
# polynomial multiplication. All multiplications use .poly_mul() instead.
.poly_from_roots <- function(roots) {
  if (length(roots) == 0L) return(1.0)
  roots <- as.complex(roots)
  n     <- length(roots)
  used  <- logical(n)
  p     <- 1.0

  for (i in seq_len(n)) {
    if (used[i]) next
    r <- roots[i]

    if (abs(Im(r)) < .Machine$double.eps^0.5) {
      # Real root: multiply by (x - Re(r))
      p       <- .poly_mul(p, c(1.0, -Re(r)))
      used[i] <- TRUE
    } else {
      # Complex root: find conjugate partner among remaining unused roots
      dists       <- abs(roots - Conj(r))
      dists[used] <- Inf
      dists[i]    <- Inf
      j           <- which.min(dists)

      if (is.finite(dists[j]) && dists[j] < 1e-8) {
        # Real quadratic: (x - r)(x - conj(r)) = x^2 - 2 Re(r) x + |r|^2
        p       <- .poly_mul(p, c(1.0, -2.0 * Re(r), Re(r)^2 + Im(r)^2))
        used[i] <- TRUE
        used[j] <- TRUE
      } else {
        # No conjugate found: complex fallback (Re() called in caller)
        p_cplx  <- c(as.complex(p), 0i) - c(0i, as.complex(p) * r)
        p       <- Re(p_cplx)
        used[i] <- TRUE
      }
    }
  }
  p
}

# Direct polynomial multiplication: p(x) * q(x).
# Coefficients stored highest-power first.
# Avoids R's convolve(), which computes cross-correlation by default (conj=TRUE).
.poly_mul <- function(p, q) {
  np  <- length(p)
  nq  <- length(q)
  res <- numeric(np + nq - 1L)
  for (i in seq_len(np)) res[i:(i + nq - 1L)] <- res[i:(i + nq - 1L)] + p[i] * q
  res
}

# scipy.signal.lfilter_zi equivalent.
# Solves for the steady-state initial conditions of a Direct Form II
# transposed filter so that both forward and backward passes of filtfilt
# start without a startup transient.
.lfilter_zi <- function(b, a) {
  b <- b / a[1L]
  a <- a / a[1L]
  n <- max(length(a), length(b))
  if (length(b) < n) b <- c(b, numeric(n - length(b)))
  if (length(a) < n) a <- c(a, numeric(n - length(a)))
  m <- n - 1L

  # Companion matrix of the denominator polynomial, then transpose
  comp <- matrix(0.0, m, m)
  comp[1L, ] <- -a[2L:n]
  if (m > 1L) {
    for (i in seq_len(m - 1L)) comp[i + 1L, i] <- 1.0
  }

  # Solve  (I - comp^T) zi = b[2:n] - a[2:n] * b[1]
  solve(diag(m) - t(comp), b[2L:n] - a[2L:n] * b[1L])
}

# Direct Form II transposed IIR filter with initial state zi.
# Equivalent to scipy.signal.lfilter(b, a, x, zi=zi).
.lfilter <- function(b, a, x, zi) {
  b <- b / a[1L]
  a <- a / a[1L]
  m <- max(length(b), length(a)) - 1L
  if (length(b) <= m) b <- c(b, numeric(m + 1L - length(b)))
  if (length(a) <= m) a <- c(a, numeric(m + 1L - length(a)))

  b1 <- b[1L]
  bt <- b[-1L]
  at <- a[-1L]

  n_x <- length(x)
  y   <- numeric(n_x)
  z   <- zi

  for (k in seq_len(n_x)) {
    yk   <- b1 * x[k] + z[1L]
    y[k] <- yk
    z_new <- bt * x[k] - at * yk
    if (m > 1L) z_new[seq_len(m - 1L)] <- z_new[seq_len(m - 1L)] + z[2L:m]
    z <- z_new
  }

  list(y = y, zf = z)
}

# Odd-reflection padding at both ends (scipy padtype = "odd").
.odd_ext <- function(x, n) {
  nx <- length(x)
  c(
    2.0 * x[1L]  - x[seq(n + 1L, 2L,      by = -1L)],
    x,
    2.0 * x[nx]  - x[seq(nx - 1L, nx - n, by = -1L)]
  )
}

# scipy.signal.filtfilt equivalent with lfilter_zi initial conditions and
# odd-reflection padding. Matches scipy output to machine precision on the
# signal body; edge agreement is also substantially improved vs zero-zi.
.filtfilt_zi <- function(b, a, x) {
  zi  <- .lfilter_zi(b, a)
  pad <- 3L * max(length(b), length(a))
  ext <- .odd_ext(x, pad)

  # Forward pass: start in steady state scaled by first padded sample
  fwd <- .lfilter(b, a, ext, zi * ext[1L])

  # Backward pass: reverse, start in steady state scaled by last fwd sample
  bwd <- .lfilter(b, a, rev(fwd$y), zi * fwd$y[length(fwd$y)])

  # Reverse and strip padding
  y <- rev(bwd$y)
  y[seq(pad + 1L, pad + length(x))]
}
