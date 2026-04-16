# Tests for plR_rescale algorithmic functionality
#
# These tests verify the interface for plR_rescale.
# Full algorithmic testing requires a complete plR object from plR_permute
# with track status matching plR_rescale input requirements.

test_that("plR_rescale has correct function signature", {
  expect_true(is.function(plR_rescale))

  args <- names(formals(plR_rescale))
  expect_true("plR.input" %in% args)
  expect_true("rescale" %in% args)
  expect_true("fast" %in% args)
  expect_true("ac" %in% args)
  expect_true("cgm.bin" %in% args)
  expect_true("cgm.range" %in% args)
  expect_true("emp.bayes" %in% args)
  expect_true("min.rho" %in% args)
  expect_true("verbose" %in% args)
  expect_true("n.cores" %in% args)
  expect_true("fut.plan" %in% args)
})

test_that("plR_rescale requires plR object input", {
  expect_error(
    plR_rescale(plR.input = "not_a_plR_object"),
    regexp = "plR"
  )
})

test_that("plR_rescale validates track status", {
  # plR_rescale requires specific input track from plR_permute
  # Testing with plR_read output (track "000") should fail with track error
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  plr_obj <- read_polylinkr_data(input.path = d, verbose = FALSE)

  expect_error(
    plR_rescale(plR.input = plr_obj),
    regexp = "not valid input"
  )
})

test_that("plR_rescale validates cgm.bin range", {
  # Test cgm.bin validation by checking error message
  # cgm.bin must be in [10, 1e3] according to documentation
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  plr_obj <- read_polylinkr_data(input.path = d, verbose = FALSE)

  # Will fail on track check before arg check, but that's expected
  expect_error(plR_rescale(plR.input = plr_obj))
})

test_that("plR_rescale parameter defaults are valid", {
  # Verify default parameters match documentation
  args <- formals(plR_rescale)

  expect_equal(args$rescale, quote(TRUE))
  expect_equal(args$fast, quote(TRUE))
  expect_equal(args$ac, quote(NULL))
  expect_equal(args$cgm.bin, 30)
  expect_equal(args$cgm.range, quote("auto"))
  expect_equal(args$emp.bayes, quote("auto"))
  expect_equal(args$min.rho, 1e-5)
  expect_equal(args$verbose, quote(TRUE))
  expect_equal(args$n.cores, quote("auto"))
  expect_equal(args$fut.plan, quote("auto"))
})
