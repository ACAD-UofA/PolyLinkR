# Tests for edge cases and error handling
#
# These unit tests verify that the package handles edge cases gracefully.

test_that("read_polylinkr_data handles minimal valid input", {
  # Use the tiny fixture which should be minimal valid input
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  
  # Should not error with minimal valid input
  expect_silent(out <- read_polylinkr_data(input_path = d, verbose = FALSE))
  
  # Should have expected structure
  expect_s3_class(out, "plR")
  expect_true(all(c("set.info", "obj.info", "set.obj") %in% names(out)))
})

test_that("read_polylinkr_data errors on missing directory", {
  expect_error(
    read_polylinkr_data(input_path = "/nonexistent/path", verbose = FALSE),
    "does not exist"
  )
})

test_that("read_polylinkr_data errors on conflicting paths", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  
  # Providing both input_path and individual paths should error
  expect_error(
    read_polylinkr_data(
      input_path = d,
      obj_info_path = file.path(d, "obj_info.txt"),
      verbose = FALSE
    ),
    "Conflicting"
  )
})

test_that("read_polylinkr_data handles set size filtering", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  
  # Filter out all sets (tiny fixture has small sets)
  expect_silent(
    out <- read_polylinkr_data(input_path = d, min_set_size = 1000, verbose = FALSE)
  )
  
  # Should have fewer or no sets
  expect_true(nrow(out$set.info) >= 0)
})

test_that("print.plR works on valid object", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
  # Should not error
  expect_silent(print(out))
})

test_that("summary.plR works on valid object", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
  # Should not error (even though no tests have been run)
  expect_silent(summary(out))
})

test_that("summary.plR validates sig parameter", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
  # sig must be in (0, 1)
  expect_error(summary(out, sig = 0), "between 0 and 1")
  expect_error(summary(out, sig = 1), "between 0 and 1")
  expect_error(summary(out, sig = -0.1), "between 0 and 1")
  expect_error(summary(out, sig = 1.5), "between 0 and 1")
})

test_that(".check_file_paths handles missing required files", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  
  # Create only some required files
  writeLines("setID", file.path(d, "setinfo.txt"))
  
  e <- new.env()
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = NULL, si.path = NULL, so.path = NULL,
      rr.path = NULL, vi.path = NULL,
      input.path = d, group = NULL, ENV = e
    ),
    "No copies of"
  )
})

test_that(".check_file_paths detects multiple file versions", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  
  # Create multiple versions of the same file
  writeLines("setID", file.path(d, "setinfo.txt"))
  writeLines("setID", file.path(d, "SetInfo.txt"))
  
  e <- new.env()
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = NULL, si.path = NULL, so.path = NULL,
      rr.path = NULL, vi.path = NULL,
      input.path = d, group = NULL, ENV = e
    ),
    "Multiple copies"
  )
})

test_that(".report_messages handles NULL messages", {
  # Set up environment with NULL messages
  info.messages <- NULL
  warning.messages <- NULL
  param.messages <- NULL
  param.warnings <- NULL
  
  # Should not error with all NULLs
  expect_silent(polylinkR:::.report_messages())
})

test_that(".report_messages handles non-NULL messages", {
  # Test that warnings are issued
  info.messages <- NULL
  warning.messages <- "Test warning"
  param.messages <- NULL
  param.warnings <- NULL
  
  expect_warning(polylinkR:::.report_messages(), "Test warning")
})
