test_that(".par_params rejects invalid fut.plan", {
  e <- new.env(parent = emptyenv())
  expect_error(
    polylinkR:::.par_params(1L, "cluster", FALSE, e),
    "fut.plan"
  )
})

test_that(".par_params sets sequential mode with one core", {
  e <- new.env(parent = emptyenv())
  polylinkR:::.par_params(1L, "sequential", FALSE, e)
  expect_equal(e$n.cores, 1L)
  expect_equal(e$fut.plan, "sequential")
})
