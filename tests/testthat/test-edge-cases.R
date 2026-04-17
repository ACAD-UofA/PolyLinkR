# Tests for edge cases and error handling
#
# These unit tests verify that the package handles edge cases gracefully.

test_that("read_polylinkr_data handles minimal valid input", {
  # Use the tiny fixture which should be minimal valid input
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  
  # Should not error with minimal valid input
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
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
      object_info_path = file.path(d, "obj_info.txt"),
      verbose = FALSE
    ),
    "Conflicting"
  )
})

test_that("read_polylinkr_data validates set size parameters", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")

  # min_set_size must be at least 2 and less than max_set_size
  expect_error(
    read_polylinkr_data(input_path = d, min_set_size = 1000, verbose = FALSE)
  )
})

test_that("print.plR produces output on valid object", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
  # Should produce output (not be silent)
  expect_output(print(out))
})

test_that("summary.plR produces output on valid object", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  out <- read_polylinkr_data(input_path = d, verbose = FALSE)
  
  # Should produce output (not be silent)
  expect_output(summary(out))
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
  
  # Create only some required files (missing objinfo and setobj)
  writeLines("setID", file.path(d, "setinfo.txt"))
  
  e <- new.env()
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = NULL, si.path = NULL, so.path = NULL,
      rr.path = NULL, vi.path = NULL,
      input.path = d, group = NULL, ENV = e
    )
  )
})

test_that(".check_file_paths detects conflicting paths", {
  d <- system.file("extdata", "tiny_polylinkr", package = "polylinkR")
  e <- new.env()
  
  # Providing both input.path and specific file paths should error
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = file.path(d, "obj_info.txt"),
      si.path = NULL, so.path = NULL,
      rr.path = NULL, vi.path = NULL,
      input.path = d, group = NULL, ENV = e
    ),
    "Conflicting"
  )
})

test_that(".create_seed returns valid seeds", {
  # NULL should generate random seed
  seed1 <- polylinkR:::.create_seed(NULL)
  expect_type(seed1, "integer")
  expect_length(seed1, 1)
  
  # Explicit seed should be returned as-is
  seed2 <- polylinkR:::.create_seed(42L)
  expect_equal(seed2, 42L)
  
  # Negative seed should work
  seed3 <- polylinkR:::.create_seed(-100L)
  expect_equal(seed3, -100L)
})

test_that(".create_seed rejects invalid seeds", {
  # Multiple values
  expect_error(polylinkR:::.create_seed(c(1, 2)), "Incorrect random seed")
  
  # Non-numeric
  expect_error(polylinkR:::.create_seed("abc"), "Incorrect random seed")
  
  # Too large
  expect_error(polylinkR:::.create_seed(.Machine$integer.max + 1), "Incorrect random seed")
})

test_that(".report_messages emits no output when all variables are NULL", {
  f <- function() {
    warning.messages <- NULL
    param.warnings   <- NULL
    info.messages    <- NULL
    param.messages   <- NULL
    withr::local_options(verbose = FALSE)
    polylinkR:::.report_messages()
  }
  expect_silent(f())
})

test_that(".report_messages issues warning for warning.messages", {
  f <- function() {
    warning.messages <- "test warning"
    param.warnings   <- NULL
    info.messages    <- NULL
    param.messages   <- NULL
    withr::local_options(verbose = FALSE)
    polylinkR:::.report_messages()
  }
  expect_warning(f(), "test warning")
})

test_that(".report_messages issues warning for param.warnings", {
  f <- function() {
    warning.messages <- NULL
    param.warnings   <- "param warning"
    info.messages    <- NULL
    param.messages   <- NULL
    withr::local_options(verbose = FALSE)
    polylinkR:::.report_messages()
  }
  expect_warning(f(), "param warning")
})

test_that(".report_messages prints info.messages when verbose", {
  f <- function() {
    warning.messages <- NULL
    param.warnings   <- NULL
    info.messages    <- "some info"
    param.messages   <- NULL
    withr::local_options(verbose = TRUE)
    polylinkR:::.report_messages()
  }
  expect_output(f(), "some info")
})

test_that(".report_messages prints param.messages when verbose", {
  f <- function() {
    warning.messages <- NULL
    param.warnings   <- NULL
    info.messages    <- NULL
    param.messages   <- "some param info"
    withr::local_options(verbose = TRUE)
    polylinkR:::.report_messages()
  }
  expect_output(f(), "some param info")
})