test_that(".create_seed returns integer for NULL user.seed", {
  seed <- polylinkR:::.create_seed(user.seed = NULL)
  expect_type(seed, "integer")
  expect_length(seed, 1L)
})

test_that(".create_seed accepts valid single integer seed", {
  expect_identical(polylinkR:::.create_seed(user.seed = 42L), 42L)
  expect_identical(polylinkR:::.create_seed(user.seed = 42), 42L)
})

test_that(".create_seed errors on invalid seed", {
  expect_error(polylinkR:::.create_seed(user.seed = "a"), "Incorrect random seed")
  expect_error(polylinkR:::.create_seed(user.seed = c(1L, 2L)), "Incorrect random seed")
  expect_error(
    polylinkR:::.create_seed(user.seed = .Machine$integer.max + 1),
    "Incorrect random seed"
  )
})

test_that(".get_time returns formatted string and difftime", {
  start <- Sys.time()
  out <- polylinkR:::.get_time(start.time = start)
  expect_type(out, "list")
  expect_length(out, 2L)
  expect_type(out[[1L]], "character")
  expect_length(out[[1L]], 1L)
  expect_s3_class(out[[2L]], "difftime")
})

test_that(".path_check errors when input.path conflicts with explicit file paths", {
  d <- withr::local_tempdir()
  e <- new.env()
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = file.path(d, "ObjInfo.txt"),
      si.path = NULL,
      so.path = NULL,
      rr.path = NULL,
      vi.path = NULL,
      input.path = d,
      group = NULL,
      ENV = e
    ),
    "Conflicting"
  )
})

test_that(".path_check errors when multi.path files are missing", {
  d <- withr::local_tempdir()
  oi <- file.path(d, "oi.txt")
  si <- file.path(d, "si.txt")
  so <- file.path(d, "so.txt")
  file.create(oi, si)
  e <- new.env()
  expect_error(
    polylinkR:::.check_file_paths(
      oi.path = oi,
      si.path = si,
      so.path = so,
      rr.path = NULL,
      vi.path = NULL,
      input.path = NULL,
      group = NULL,
      ENV = e
    ),
    "does not exist"
  )
})

test_that(".plR_track returns expected data.table columns", {
  tab <- polylinkR:::.plR_track()
  expect_s3_class(tab, "data.table")
  expect_named(tab, c("FUNCTION", "INPUT", "OUTPUT"))
  expect_identical(
    tab[["FUNCTION"]],
    c("read_polylinkr_data", "permute_polylinkr_data", "rescale_polylinkr_data", "prune_polylinkr_data")
  )
})
