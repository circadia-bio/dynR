# Suppress "no visible binding" R CMD check notes for ggplot2 aesthetics.
utils::globalVariables(c("row", "col", "val", "timepoint", "R", "state", "rss"))

# Diverging dynR palette for FC matrices: deep indigo -> periwinkle -> brick red.
# Matches the corr_corr vignette: colorRampPalette(c("#341C5D","#E8ECF8","#9E3C30")).
.fc_palette <- c("#341C5D", "#8D9FD7", "#E8ECF8", "#E19A8F", "#9E3C30")

#' Plot a functional connectivity matrix
#'
#' Render a square correlation or phase-locking matrix as a heatmap using the
#' dynR diverging palette (deep indigo -> periwinkle -> brick red), centred at
#' zero. Matches the colour scheme used in the dynR vignettes and pkgdown site.
#'
#' @param fc_matrix Numeric matrix \[N x N\]. Pearson correlations or
#'   phase-locking values.
#' @param title Character. Plot title. Default `"Functional connectivity"`.
#' @param limits Numeric vector `c(lo, hi)`. Colour scale limits.
#'   Default `c(-1, 1)`.
#'
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_raster scale_fill_gradientn coord_fixed labs theme_minimal theme element_blank element_text
#' @export
#' @examples
#' data(fc, package = "dynR")
#' plot_fc(fc[1:20, 1:20])
plot_fc <- function(fc_matrix, title = "Functional connectivity",
                    limits = c(-1, 1)) {
  N  <- nrow(fc_matrix)
  df <- data.frame(
    row = rep(seq_len(N), N),
    col = rep(seq_len(N), each = N),
    val = as.vector(fc_matrix)
  )
  ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = val)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colours = .fc_palette,
      limits  = limits,
      name    = "r"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "Parcel", y = "Parcel", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(face = "bold")
    )
}

#' Plot Kuramoto synchrony time series
#'
#' Line plot of the Kuramoto order parameter R(t) over time, with a dashed
#' horizontal line at the mean. Uses the dynR periwinkle/brick-red palette.
#'
#' @param synchrony Numeric vector. Kuramoto R(t), as returned by `kuramoto()`
#'   or the `$synchrony` element of `leida_pipeline()`.
#' @param title Character. Plot title.
#'
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_y_continuous labs theme_minimal theme element_blank element_text
#' @importFrom stats sd
#' @export
#' @examples
#' set.seed(1)
#' ts <- matrix(rnorm(10 * 200), nrow = 10)
#' ph  <- hilbert_phases(ts)
#' kop <- kuramoto(ph)
#' plot_synchrony(kop$synchrony)
plot_synchrony <- function(synchrony, title = "Kuramoto order parameter") {
  df <- data.frame(timepoint = seq_along(synchrony), R = synchrony)
  mn <- mean(synchrony)
  ggplot2::ggplot(df, ggplot2::aes(x = timepoint, y = R)) +
    ggplot2::geom_line(colour = "#8D9FD7", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = mn, colour = "#9E3C30",
                        linetype = "dashed", linewidth = 0.9) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x        = "Timepoint",
      y        = "R(t)",
      title    = title,
      subtitle = paste0(
        "Metastability = ", round(stats::sd(synchrony), 4),
        "   |   Mean = ",   round(mn, 4)
      )
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold")
    )
}

#' Plot a brain state sequence
#'
#' Tile plot of brain state labels over time (or windows), coloured with the
#' dynR state palette (periwinkle, brick red, mauve, deep indigo, dusty rose).
#'
#' @param states Integer or factor vector of state labels.
#' @param x_label Character. x-axis label. Default `"Timepoint"`.
#' @param title Character. Plot title. Default `"Brain state sequence"`.
#' @param palette Character vector of colours. Defaults to the dynR 5-colour
#'   state palette, recycled if `K > 5`.
#'
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual labs theme_minimal theme element_blank element_text
#' @export
#' @examples
#' set.seed(1)
#' states <- sample(1:4, 100, replace = TRUE)
#' plot_state_sequence(states)
plot_state_sequence <- function(states,
                                x_label = "Timepoint",
                                title   = "Brain state sequence",
                                palette = NULL) {
  K   <- length(unique(states))
  pal <- if (is.null(palette)) .state_cols(K) else rep_len(palette, K)
  df  <- data.frame(timepoint = seq_along(states), state = factor(states))
  ggplot2::ggplot(df, ggplot2::aes(x = timepoint, y = 1, fill = state)) +
    ggplot2::geom_tile(height = 1) +
    ggplot2::scale_fill_manual(values = pal, name = "State") +
    ggplot2::labs(x = x_label, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(face = "bold")
    )
}

#' Plot method for dynR_leida objects
#'
#' @param x A `dynR_leida` object as returned by `leida_pipeline()`.
#' @param type Character. `"synchrony"` (default) plots the Kuramoto order
#'   parameter time series via `plot_synchrony()`; `"fc"` plots the
#'   time-averaged phase-locking matrix via `plot_fc()`.
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A `ggplot` object.
#' @export
plot.dynR_leida <- function(x, type = c("synchrony", "fc"), ...) {
  type <- match.arg(type)
  if (type == "synchrony") {
    plot_synchrony(x$synchrony, ...)
  } else {
    mean_fc <- apply(x$sync_conn, c(1L, 2L), mean)
    plot_fc(mean_fc, title = "Mean phase-locking matrix", ...)
  }
}

#' Plot method for dynR_sw objects
#'
#' @param x A `dynR_sw` object as returned by `sw_pipeline()`.
#' @param type Character. `"rss"` (default) plots the RSS cofluctuation time
#'   series; `"fc"` plots the time-averaged sliding-window FC matrix via
#'   `plot_fc()`.
#' @param ... Additional arguments passed to the underlying plot function.
#'
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot aes geom_line geom_hline labs theme_minimal theme element_blank element_text
#' @importFrom stats sd
#' @export
plot.dynR_sw <- function(x, type = c("rss", "fc"), ...) {
  type <- match.arg(type)
  if (type == "rss") {
    thr <- mean(x$rss) + 2 * stats::sd(x$rss)
    df  <- data.frame(timepoint = seq_along(x$rss), rss = x$rss)
    ggplot2::ggplot(df, ggplot2::aes(x = timepoint, y = rss)) +
      ggplot2::geom_line(colour = "#8D9FD7", linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = thr, colour = "#9E3C30",
                          linetype = "dashed", linewidth = 0.9) +
      ggplot2::labs(
        x        = "Timepoint",
        y        = "RSS",
        title    = "Root-sum-square cofluctuation",
        subtitle = paste0(round(mean(x$rss > thr) * 100L, 1L),
                          "% of timepoints above mean + 2 SD")
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(face = "bold")
      )
  } else {
    mean_fc <- apply(x$corr_mats, c(1L, 2L), mean)
    plot_fc(mean_fc, title = "Mean sliding-window FC", ...)
  }
}
