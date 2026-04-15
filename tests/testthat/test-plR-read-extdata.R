test_that("plR_read succeeds on bundled tiny_polylinkR fixture", {
  d <- system.file("extdata", "tiny_polylinkR", package = "polylinkR")
  skip_if(d == "", "tiny_polylinkR fixture not installed")

  out <- polylinkR::plR_read(
    input.path = d,
    verbose = FALSE,
    min.set.n = 2L,
    max.set.n = 100L,
    set.merge = 0.95
  )
  expect_s3_class(out, "plR")
  expect_true(nrow(out$obj.info) >= 2L)
  expect_true(nrow(out$set.info) >= 1L)
})
