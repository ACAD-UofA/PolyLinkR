test_that(".path_check errors when input.path is not a directory", {
  d <- withr::local_tempdir()
  e <- new.env()
  expect_error(
    polylinkR:::.path_check(
      oi.path = NULL,
      si.path = NULL,
      so.path = NULL,
      rr.path = NULL,
      vi.path = NULL,
      input.path = file.path(d, "missing-dir"),
      group = NULL,
      ENV = e
    ),
    "directory does not exist"
  )
})

test_that(".path_check errors when separate paths are incomplete", {
  e <- new.env()
  expect_error(
    polylinkR:::.path_check(
      oi.path = "/nonexistent/obj.txt",
      si.path = NULL,
      so.path = NULL,
      rr.path = NULL,
      vi.path = NULL,
      input.path = NULL,
      group = NULL,
      ENV = e
    ),
    "Incomplete"
  )
})
