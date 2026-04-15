# Integration test: keep local and CRAN fast; skip on CI agents.
test_that("plR_read succeeds on a small HumanPops extract", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()

  oi_full <- system.file("data", "HumanPops_ObjInfo.txt", package = "polylinkR")
  si_full <- system.file("data", "HumanPops_SetInfo.txt", package = "polylinkR")
  so_full <- system.file("data", "HumanPops_SetObj.txt", package = "polylinkR")
  skip_if(oi_full == "", "Bundled HumanPops extract missing")

  d <- withr::local_tempdir()
  oi <- utils::read.delim(oi_full, nrows = 120L, sep = "\t", stringsAsFactors = FALSE)
  utils::write.table(oi, file.path(d, "ObjInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  so <- utils::read.delim(so_full, sep = "\t", stringsAsFactors = FALSE)
  so <- so[so$objID %in% oi$objID, , drop = FALSE]
  so <- so[so$setID == min(so$setID), , drop = FALSE]
  utils::write.table(so, file.path(d, "SetObj.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  si <- utils::read.delim(si_full, sep = "\t", stringsAsFactors = FALSE)
  si <- si[si$setID %in% unique(so$setID), , drop = FALSE]
  utils::write.table(si, file.path(d, "SetInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  out <- polylinkR::plR_read(
    input.path = d,
    verbose = FALSE,
    min.set.n = 2L,
    max.set.n = Inf,
    set.merge = 0.95
  )
  expect_s3_class(out, "plR")
  expect_true(all(c("set.info", "obj.info", "set.obj") %in% names(out)))
})
