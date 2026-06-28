# data-raw/generate_test_fixtures.R
#
# Run once to generate RDS test fixtures from the dynfc Python .npy reference
# data. Requires the dynfc repo to be checked out alongside dynR at
# ../dynfc/data/.
#
# Usage (from dynR project root):
#   source("data-raw/generate_test_fixtures.R")
#
# Requires: reticulate, numpy (in the active Python env)

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
ts <- np$load(file.path(dynfc_data, "ts.npy"))   # [200, 600]
fc <- np$load(file.path(dynfc_data, "fc.npy"))   # [200, 200]

# 10-parcel subset (rows 1–10 in R = 0–9 in Python)
ts10 <- ts[1:10, ]     # [10, 600]
fc10 <- fc[1:10, 1:10] # [10, 10]

# Compute ets10: z-score each parcel row (ddof = 1), then element-wise
# product for every unique parcel pair (upper triangle, k = 1)
ts10_z <- t(scale(t(ts10)))
idx    <- which(upper.tri(matrix(0, 10, 10)), arr.ind = TRUE)  # 45 edges
ets10  <- ts10_z[idx[, 1], ] * ts10_z[idx[, 2], ]             # [45, 600]

message("Saving fixtures to inst/testdata/...")
dir.create("inst/testdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(ts10,  "inst/testdata/ts10.rds")
saveRDS(fc10,  "inst/testdata/fc10.rds")
saveRDS(ets10, "inst/testdata/ets10.rds")
message("Done. Fixtures written:")
message("  inst/testdata/ts10.rds  — [10 x 600] BOLD timeseries (10 parcels)")
message("  inst/testdata/fc10.rds  — [10 x 10] FC matrix (ground truth)")
message("  inst/testdata/ets10.rds — [45 x 600] edge time series (ground truth)")
