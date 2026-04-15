minimal_plr <- function(track = "000") {
  SI <- data.table::data.table(setID = 1L, setSize = 2L)
  OI <- data.table::data.table(objID = 1:2, objStat = c(0.1, 0.2))
  SO <- data.table::data.table(setID = 1L, objID = 1:2)
  BASE <- list(set.info = SI, obj.info = OI, set.obj = SO)
  x <- polylinkR:::.new_plR(BASE = BASE, plR.args = NULL, plR.data = NULL)
  attr(x, "plR.track") <- track
  x
}

test_that("print.plR reports empty object", {
  empty <- polylinkR:::.new_plR(BASE = NA)
  expect_output(print(empty), "Empty plR object")
})

test_that("summary.plR rejects sig outside (0, 1)", {
  x <- minimal_plr()
  expect_error(summary(x, sig = 0), "between 0 and 1")
  expect_error(summary(x, sig = 1), "between 0 and 1")
})

test_that("summary.plR reports when no test columns exist", {
  x <- minimal_plr()
  expect_output(summary(x, sig = 0.05), "No tests performed")
})

test_that(".plR_check errors for non-plR object", {
  f <- function(plR.input) {
    polylinkR:::.plR_check("permute", environment())
  }
  expect_error(f(list(a = 1)), "not a plR class object")
})

test_that(".plR_check errors when track is invalid for target function", {
  f <- function(plR.input) {
    polylinkR:::.plR_check("permute", environment())
  }
  x <- minimal_plr("999")
  expect_error(f(x), "not valid input for plR_permute")
})

test_that(".plR_check succeeds for permute with initial read track", {
  f <- function(plR.input) {
    polylinkR:::.plR_check("permute", environment())
    plr.track
  }
  x <- minimal_plr("000")
  expect_equal(unname(f(x)), c(0, 0, 0), tolerance = 0)
})
