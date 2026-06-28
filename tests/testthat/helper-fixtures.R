# Helper: load a test fixture from inst/testdata/, skipping if absent.
# Used by all test files that require reference data from dynfc.
load_fixture <- function(name) {
  path <- system.file("testdata", name, package = "dynR")
  if (nchar(path) == 0 || !file.exists(path)) {
    testthat::skip(paste0(
      "Test fixture '", name, "' not found. ",
      "Run data-raw/generate_test_fixtures.R to generate it."
    ))
  }
  readRDS(path)
}
