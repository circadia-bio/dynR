# data-raw/generate_test_fixtures.R
#
# Run once to generate:
#   - data/ts.rda, data/fc.rda       — full 200-parcel package data (vignettes, examples)
#   - inst/testdata/ts10.rds          — 10-parcel BOLD subset (tests)
#   - inst/testdata/fc10.rds          — 10-parcel FC matrix (tests)
#   - inst/testdata/ets10.rds         — 10-parcel edge time series (tests)
#
# Requires the dynfc repo at ../dynfc/data/ and reticulate + numpy.
#
# Usage (from dynR project root):
#   source("data-raw/generate_test_fixtures.R")

library(reticulate)

np <- import("numpy")

dynfc_data <- "../dynfc/data"
if (!dir.exists(dynfc_data)) {
  stop(
    "dynfc data directory not found at ", dynfc_data, "\n",
    "Clone https://github.com/LucasFranca/dynfc alongside dynR."
  )
}

message("Loading .npy reference data...")
ts <- np$load(file.path(dynfc_data, "ts.npy"))  # [200, 600]
fc <- np$load(file.path(dynfc_data, "fc.npy"))  # [200, 200]

# ── Full package data ────────────────────────────────────────────────────────
message("Saving full package data to data/...")
dir.create("data", showWarnings = FALSE)
save(ts, file = "data/ts.rda", compress = "bzip2")
save(fc, file = "data/fc.rda", compress = "bzip2")
message("  data/ts.rda — [200 x 600] BOLD timeseries (200 parcels, 600 timepoints)")
message("  data/fc.rda — [200 x 200] functional connectivity matrix")

# ── 10-parcel test fixtures ──────────────────────────────────────────────────
message("Generating 10-parcel test fixtures...")
ts10 <- ts[1:10, ]
fc10 <- fc[1:10, 1:10]

ts10_z <- t(scale(t(ts10)))
idx    <- which(upper.tri(matrix(0, 10, 10)), arr.ind = TRUE)  # 45 edges
ets10  <- ts10_z[idx[, 1], ] * ts10_z[idx[, 2], ]             # [45, 600]

message("Saving test fixtures to inst/testdata/...")
dir.create("inst/testdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(ts10,  "inst/testdata/ts10.rds")
saveRDS(fc10,  "inst/testdata/fc10.rds")
saveRDS(ets10, "inst/testdata/ets10.rds")
message("  inst/testdata/ts10.rds  — [10 x 600] BOLD timeseries (10 parcels)")
message("  inst/testdata/fc10.rds  — [10 x 10] FC matrix")
message("  inst/testdata/ets10.rds — [45 x 600] edge time series")

message("\nDone.")
