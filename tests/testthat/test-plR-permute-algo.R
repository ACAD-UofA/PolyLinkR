# Tests for plR_permute function interface

test_that("plR_permute has correct function signature", {
  expect_true(is.function(plR_permute))

  # Verify key arguments exist
  args <- names(formals(plR_permute))
  expect_true("plR.input" %in% args)
  expect_true("permute" %in% args)
  expect_true("n.perm" %in% args)
  expect_true("n.boot" %in% args)
  expect_true("alt" %in% args)
  expect_true("gpd.cutoff" %in% args)
  expect_true("seed" %in% args)
})

test_that("plR_permute validates alt parameter", {
  # Use tiny fixture for valid plR object
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  plr_obj <- read_polylinkr_data(input_path = d, verbose = FALSE)

  # Invalid alt should error
  expect_error(
    plR_permute(plR.input = plr_obj, alt = "invalid"),
    regexp = "alt"
  )
})

test_that("plR_permute requires plR object input", {
  expect_error(
    plR_permute(plR.input = "not_a_plR_object"),
    regexp = "plR"
  )
})

test_that("plR_permute n.perm validation", {
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  plr_obj <- read_polylinkr_data(input.path = d, verbose = FALSE)

  # n.perm must be >= 1e5 and divisible by 1e5
  expect_error(
    plR_permute(plR.input = plr_obj, n.perm = 1000),
    regexp = "n.perm"
  )
})
