# Tests for permute_polylinkr_data function interface

test_that("permute_polylinkr_data has correct function signature", {
  expect_true(is.function(permute_polylinkr_data))

  # Verify key arguments exist
  args <- names(formals(permute_polylinkr_data))
  expect_true("plr_input" %in% args)
  expect_true("permute" %in% args)
  expect_true("n_permutations" %in% args)
  expect_true("n_bootstraps" %in% args)
  expect_true("alternative" %in% args)
  expect_true("gpd_cutoff" %in% args)
  expect_true("seed" %in% args)
})

test_that("permute_polylinkr_data validates alt parameter", {
  # Use tiny fixture for valid polylinkr object
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  polylinkr_obj <- read_polylinkr_data(input_path = d, verbose = FALSE)

  # Invalid alt should error
  expect_error(
    permute_polylinkr_data(plr_input = polylinkr_obj, alt = "invalid"),
    regexp = "alt"
  )
})

test_that("permute_polylinkr_data requires polylinkr object input", {
  expect_error(
    permute_polylinkr_data(plr_input = "not_a_polylinkr_object"),
    regexp = "plr"
  )
})

test_that("permute_polylinkr_data n.perm validation", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  polylinkr_obj <- read_polylinkr_data(input_path = d, verbose = FALSE)

  # n.perm must be >= 1e5 and divisible by 1e5
  expect_error(
    permute_polylinkr_data(plr_input = polylinkr_obj, n.perm = 1000),
    regexp = "n.perm"
  )
})
