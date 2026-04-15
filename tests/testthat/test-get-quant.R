test_that(".get_quant returns aligned quantiles and error band", {
  skip_if_not_installed("tdigest")
  d1 <- tdigest::td_create()
  d2 <- tdigest::td_create()
  for (x in stats::runif(80L)) {
    d1 <- tdigest::td_add(d1, x, 1L)
    d2 <- tdigest::td_add(d2, x + 0.01, 1L)
  }
  tX <- list(d1, d2)
  out <- polylinkR:::.get_quant(tX, eQnt = c(0.25, 0.5), eps = 1e-4)
  expect_type(out, "list")
  expect_named(out, c("qX", "qE"))
  expect_length(out$qX, 2L)
  expect_true(all(is.finite(out$qX)))
  expect_true(all(is.finite(out$qE)))
})
