test_that(".md2p returns row-normalised weights for Gaussian kernel", {
  wt <- matrix(c(0.1, 0.5, 1, 2), nrow = 2, ncol = 2)
  p <- polylinkR:::.md2p(wt, kern.scale = 1, kern.func = "gaussian")
  expect_equal(rowSums(p), c(1, 1), tolerance = 1e-10)
})

test_that(".md2p returns row-normalised weights for exponential kernel", {
  wt <- matrix(c(0.1, 0.5, 1, 2), nrow = 2, ncol = 2)
  p <- polylinkR:::.md2p(wt, kern.scale = 1, kern.func = "exponential")
  expect_equal(rowSums(p), c(1, 1), tolerance = 1e-10)
  expect_true(all(p > 0))
})

test_that(".md2p handles inverse kernel", {
  wt <- matrix(1:4, nrow = 2, ncol = 2)
  p <- polylinkR:::.md2p(wt, kern.scale = 2, kern.func = "inverse")
  expect_equal(rowSums(p), c(1, 1), tolerance = 1e-10)
})

test_that(".cap_probs respects maxP and row sums", {
  mm <- matrix(0.5, nrow = 3, ncol = 4)
  mm <- mm / rowSums(mm)
  out <- polylinkR:::.cap_probs(mm, maxP = 0.2)
  expect_true(all(out <= 0.2 + 1e-10))
  expect_equal(rowSums(out), rep(1, 3), tolerance = 1e-10)
})

test_that(".cov_fun returns positive finite values for interior distances", {
  h <- c(0.1, 1, 5)
  v <- polylinkR:::.cov_fun(h, ovlp = 0.5, beta1 = 0.2, beta2 = 0.1, scale = 2)
  expect_true(all(is.finite(v)))
  expect_true(all(v > 0))
})

test_that(".facet_opt returns sensible layout for small panels", {
  o <- polylinkR:::.facet_opt(n.panel = 4L, max.facets = 12L)
  expect_named(o, c("n.row", "n.col", "n.pg"))
  expect_true(o$n.row * o$n.col <= 12L)
})

test_that(".args2title drops nulls and formats", {
  t <- polylinkR:::.args2title(list(a = 1, b = 2, seed = 3L, input.path = "x"))
  expect_type(t, "character")
  expect_length(t, 1L)
  expect_match(t, "Parameters")
})

test_that(".fit_gpd returns NA vector on constant input", {
  out <- polylinkR:::.fit_gpd(rep(1, 20))
  expect_length(out, 4L)
  expect_true(all(is.na(out)))
})
