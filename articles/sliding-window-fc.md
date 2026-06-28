# Correlation-based dynamic FC: sliding windows and edge cofluctuations

## Overview

Correlation-based dynFC methods characterise how functional connectivity
evolves over time by either **windowing** the timeseries into short
segments and computing Pearson correlations within each window, or by
analysing the **instantaneous co-activation** of parcel pairs as a
continuous edge time series.

`dynR` implements three correlation-based functions:

| Function | Output | Method |
|----|----|----|
| [`corr_slide()`](https://CoDe-Neuro.github.io/dynR/reference/corr_slide.md) | FC matrices 
``` math
N × N × windows
``` | Sliding-window Pearson correlation |
| [`cofluct()`](https://CoDe-Neuro.github.io/dynR/reference/cofluct.md) | Edge time series 
``` math
n\_edges × Tmax
```
 + RSS | Edge-centric cofluctuation |
| [`corr_corr()`](https://CoDe-Neuro.github.io/dynR/reference/corr_corr.md) | Correlation-of-correlations 
``` math
Tmax × Tmax
``` | Hansen et al. (2015) |

------------------------------------------------------------------------

## Sliding-window correlation

[`corr_slide()`](https://CoDe-Neuro.github.io/dynR/reference/corr_slide.md)
divides the timeseries into overlapping or non-overlapping windows and
computes the Pearson FC matrix within each one. The choice of window
length is a key parameter: too short and the estimate is noisy; too long
and rapid transitions are blurred.

``` r

# 30-timepoint windows, step of 5 (overlapping)
sw <- corr_slide(ts, window = 30, step = 5)

cat("Window size: 30 timepoints\n")
#> Window size: 30 timepoints
cat("Step:        5 timepoints\n")
#> Step:        5 timepoints
cat("N windows:  ", dim(sw$corr_mats)[3], "\n")
#> N windows:   115
```

### Validating against static FC

When `window = ncol(ts)` the single window is the full timeseries, so
the result must equal the static FC matrix.

``` r

full <- corr_slide(ts, window = ncol(ts))
max_diff <- max(abs(full$corr_mats[, , 1] - fc))
cat("Max deviation from static FC:", formatC(max_diff, format = "e"), "\n")
#> Max deviation from static FC: 8.8818e-16
```

### FC variability across windows

A useful summary is how much the inter-parcel correlations vary across
windows. High variability means the connectivity structure is genuinely
non-stationary.

``` r

# Upper triangle indices (excluding diagonal)
n_parcels <- dim(sw$corr_mats)[1]
ut <- which(upper.tri(diag(n_parcels)), arr.ind = TRUE)

# Extract upper triangle for each window → [n_edges x n_windows] matrix
edge_by_window <- apply(sw$corr_mats, 3, function(m) m[ut])

# Mean standard deviation across edges
mean_sd <- mean(apply(edge_by_window, 1, sd))
cat("Mean edge SD across windows:", round(mean_sd, 4), "\n")
#> Mean edge SD across windows: 0.5091
```

------------------------------------------------------------------------

## Edge-centric cofluctuations

The edge-centric framework (Esfahlani et al., 2020; Faskowitz et al.,
2020) does not window the data at all. Instead, it z-scores each
parcel’s timeseries and computes the element-wise product for every
unique parcel pair:

``` math
\text{ets}_{ij}(t) = z_i(t) \cdot z_j(t)
```

The result is an **edge time series** — a continuous, frame-by-frame
estimate of co-activation for every pair.

``` r

ec <- cofluct(ts)

n_edges <- nrow(ec$edge_ts)
cat("Edges (N*(N-1)/2):", n_edges, "\n")   # 200*199/2 = 19900
#> Edges (N*(N-1)/2): 19900
cat("Timepoints:       ", ncol(ec$edge_ts), "\n")
#> Timepoints:        600
```

### Root-sum-square (RSS) cofluctuation

The RSS vector summarises the total cofluctuation amplitude at each
timepoint:

``` math
\text{RSS}(t) = \sqrt{\sum_{(i,j)} \text{ets}_{ij}(t)^2}
```

High-RSS frames are moments of unusually strong co-activation across the
brain and have been shown to drive the structure of the static FC matrix
(Esfahlani et al., 2020).

``` r

plot(ec$rss, type = "l", col = "#1B6799",
     xlab = "Timepoint", ylab = "RSS",
     main = "Root-sum-square cofluctuation",
     lwd = 1.2)
abline(h = mean(ec$rss) + 2 * sd(ec$rss),
       col = "#FC544A", lty = 2, lwd = 1.5)
legend("topright", legend = "Mean + 2 SD",
       col = "#FC544A", lty = 2, lwd = 1.5, bty = "n")
```

![](sliding-window-fc_files/figure-html/rss-plot-1.png)

### High-amplitude events

Frames with RSS above mean + 2 SD are candidate **high-amplitude
cofluctuation events** — the moments that most strongly shape static FC.

``` r

threshold  <- mean(ec$rss) + 2 * sd(ec$rss)
high_frames <- which(ec$rss > threshold)
cat("High-amplitude frames (RSS > mean + 2SD):", length(high_frames), "\n")
#> High-amplitude frames (RSS > mean + 2SD): 20
cat("Proportion of scan:", round(length(high_frames) / ncol(ts), 3), "\n")
#> Proportion of scan: 0.033
```

------------------------------------------------------------------------

## Correlation of correlations

[`corr_corr()`](https://CoDe-Neuro.github.io/dynR/reference/corr_corr.md)
computes the full
``` math
Tmax × Tmax
```
matrix of pairwise correlations between edge time series (Hansen et al.,
2015). Entry *(t1, t2)* reflects how similar the co-activation patterns
were at timepoints *t1* and *t2* — a non-windowed measure of FC
recurrence.

``` r

cc <- corr_corr(ts)
dim(cc)   # [600 x 600]
#> [1] 600 600
```

The matrix is symmetric by construction and has ones on the diagonal.

``` r

cat("Symmetric:   ", isTRUE(all.equal(cc, t(cc))), "\n")
#> Symmetric:    TRUE
cat("Diagonal == 1:", all(abs(diag(cc) - 1) < 1e-10), "\n")
#> Diagonal == 1: TRUE
```

------------------------------------------------------------------------

## Preparing for brain state analysis

Both sliding-window FC and the correlation-of-correlations matrix can be
used for K-means brain state discovery. The standard approach with
sliding-window FC is to vectorise the upper triangle of each window’s FC
matrix and cluster the resulting feature vectors:

``` r

# Vectorise upper triangle for each window
ut       <- which(upper.tri(diag(n_parcels)), arr.ind = TRUE)
features <- t(apply(sw$corr_mats, 3, function(m) m[ut]))
# features: [n_windows x n_edges]

set.seed(42)
K  <- 5
km <- kmeans(features, centers = K, nstart = 50, iter.max = 500)
table(km$cluster)
#> 
#>  1  2  3  4  5 
#> 16 50 14 11 24
```

The state sequence can then be passed to `stateR` for fractional
occupancy, dwell time, and Markov transition analysis.

------------------------------------------------------------------------

## References

Hansen, E. C. A. et al. (2015). Functional connectivity dynamics:
Modeling the switching behavior of the resting state. *NeuroImage*, 105,
525–535. <https://doi.org/10.1016/j.neuroimage.2014.11.001>

Esfahlani, F. Z. et al. (2020). High-amplitude cofluctuations in
cortical activity drive functional connectivity. *PNAS*, 117(45),
28393–28401. <https://doi.org/10.1073/pnas.2005531117>

Faskowitz, J. et al. (2020). Edge-centric functional network
representations of human cerebral cortex reveal overlapping system-level
architecture. *Nature Neuroscience*, 23(12), 1644–1654.
<https://doi.org/10.1038/s41593-020-00719-y>
