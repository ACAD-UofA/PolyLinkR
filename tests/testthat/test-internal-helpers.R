# Tests for internal helper functions used in plR_permute, plR_rescale, plR_prune
#
# These unit tests verify the core statistical helpers work correctly.

test_that(".get_time returns formatted time string", {
  start_time <- Sys.time()
  Sys.sleep(0.1)
  result <- polylinkR:::.get_time(start_time)

  expect_type(result, "list")
  expect_length(result, 2)
  expect_type(result[[1]], "character")
  expect_s3_class(result[[2]], "difftime")
})

test_that(".plR_track returns valid track table", {
  tab <- polylinkR:::.plR_track()

  expect_s3_class(tab, "data.table")
  expect_named(tab, c("FUNCTION", "INPUT", "OUTPUT"))
  expect_equal(tab$FUNCTION,
    c("plR_read", "plR_permute", "plR_rescale", "plR_prune"))
})

test_that(".plR_track defines valid input/output transitions", {
  tab <- polylinkR:::.plR_track()

  # plR_read outputs "000"
  expect_true(grepl("000", tab$OUTPUT[1]))

  # plR_permute accepts "000"
  expect_true(grepl("000", tab$INPUT[2]))

  # Each function has defined output (at least one value)
  expect_true(all(nzchar(tab$OUTPUT)))
})

test_that("plR_read returns valid plR object structure", {
  # Use tiny fixture
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  out <- plR_read(input.path = d, verbose = FALSE)

  expect_s3_class(out, "plR")
  expect_true("set.info" %in% names(out))
  expect_true("obj.info" %in% names(out))
  expect_true("set.obj" %in% names(out))

  # Check required attributes
  expect_true("plR.track" %in% names(attributes(out)))
  expect_equal(attr(out, "plR.track"), "000")
})

test_that("plR_read validates min.set.n and max.set.n", {
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")

  # max.set.n too small
  expect_error(
    plR_read(input.path = d, max.set.n = 1),
    regexp = "max.set.n"
  )
})

test_that("plR_read set.merge parameter validation", {
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")

  # set.merge must be in (0, 1]
  expect_error(
    plR_read(input.path = d, set.merge = 0),
    regexp = "set.merge"
  )
  expect_error(
    plR_read(input.path = d, set.merge = 1.1),
    regexp = "set.merge"
  )
})
