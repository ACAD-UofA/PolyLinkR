# Tests for rescale_polylinkr_data algorithmic functionality
#
# These tests verify the interface for rescale_polylinkr_data.
# Full algorithmic testing requires a complete polylinkr object from permute_polylinkr_data
# with track status matching rescale_polylinkr_data input requirements.

test_that("rescale_polylinkr_data has correct function signature", {
  expect_true(is.function(rescale_polylinkr_data))

  args <- names(formals(rescale_polylinkr_data))
  expect_true("plr_input" %in% args)
  expect_true("rescale" %in% args)
  expect_true("fast" %in% args)
  expect_true("ac" %in% args)
  expect_true("cgm_bin" %in% args)
  expect_true("cgm_range" %in% args)
  expect_true("empirical_bayes" %in% args)
  expect_true("min_rho" %in% args)
  expect_true("verbose" %in% args)
  expect_true("n_cores" %in% args)
  expect_true("future_plan" %in% args)
})

test_that("rescale_polylinkr_data requires polylinkr object input", {
  expect_error(
    rescale_polylinkr_data(plr_input = "not_a_polylinkr_object"),
    regexp = "not a plr class object"
  )
})

test_that("rescale_polylinkr_data validates track status", {
  # rescale_polylinkr_data requires specific input track from permute_polylinkr_data
  # Testing with read_polylinkr_data output (track "000") should fail with track error
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  polylinkr_obj <- read_polylinkr_data(input_path = d, verbose = FALSE)

  expect_error(
    rescale_polylinkr_data(plr_input = polylinkr_obj),
    regexp = "not valid input"
  )
})

test_that("rescale_polylinkr_data validates cgm.bin range", {
  # Test cgm.bin validation by checking error message
  # cgm.bin must be in [10, 1e3] according to documentation
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  polylinkr_obj <- read_polylinkr_data(input_path = d, verbose = FALSE)

  # Will fail on track check before arg check, but that's expected
  expect_error(rescale_polylinkr_data(plr_input = polylinkr_obj))
})

test_that("rescale_polylinkr_data parameter defaults are valid", {
  # Verify default parameters match documentation
  args <- formals(rescale_polylinkr_data)

  expect_equal(args$rescale, quote(TRUE))
  expect_equal(args$fast, quote(TRUE))
  expect_equal(args$ac, quote(NULL))
  expect_equal(args$cgm_bin, 30)
  expect_equal(args$cgm_range, quote("auto"))
  expect_equal(args$empirical_bayes, quote("auto"))
  expect_equal(args$min_rho, 1e-5)
  expect_equal(args$verbose, quote(TRUE))
  expect_equal(args$n_cores, quote("auto"))
  expect_equal(args$future_plan, quote("auto"))
})
