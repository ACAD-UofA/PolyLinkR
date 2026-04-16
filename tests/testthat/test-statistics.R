# Tests for statistical functions
#
# These unit tests verify the core statistical helpers work correctly.

test_that(".fit_gpd_tail handles valid data", {
  # Generate some valid data from an exponential distribution
  set.seed(42)
  x <- rexp(1000, rate = 2)
  
  result <- polylinkR:::.fit_gpd_tail(x)
  
  # Should return 4 values: mle (2) + se^2 (2)
  expect_length(result, 4)
  
  # Parameters should be finite (not NA if fitting succeeded)
  expect_true(all(is.finite(result)))
})

test_that(".fit_gpd_tail handles edge cases", {
  # Single value - should return NAs
  result <- polylinkR:::.fit_gpd_tail(c(1))
  expect_length(result, 4)
  expect_true(all(is.na(result)))
})

test_that(".get_pvalue_upper returns valid p-values", {
  # Create mock ECDF and GPD objects
  set.seed(42)
  null_scores <- rnorm(1000)
  ecdf_fn <- ecdf(null_scores)
  
  # Create mock GPD parameters
  gpd_params <- data.table::data.table(
    shape = 0.1,
    scale = 0.5,
    exc.z = 1.5,
    cut.z = 2.0,
    min.z = 3.0,
    min.p = 1e-4,
    adj.p = 1.0
  )
  
  # Test p-values for scores in the ECDF range
  test_scores <- c(0, 0.5, 1.0, 1.5)
  pvals <- polylinkR:::.get_pvalue_upper(test_scores, 10, gpd_params, list(ecdf_fn))
  
  # P-values should be in [0, 1]
  expect_true(all(pvals >= 0 & pvals <= 1))
  
  # P-values should be monotonically decreasing with score
  expect_true(all(diff(pvals) <= 0))
})

test_that(".get_pvalue_lower returns valid p-values", {
  # Create mock ECDF and GPD objects
  set.seed(42)
  null_scores <- rnorm(1000)
  ecdf_fn <- ecdf(null_scores)
  
  # Create mock GPD parameters (negated for lower tail)
  gpd_params <- data.table::data.table(
    shape = 0.1,
    scale = 0.5,
    exc.z = -1.5,
    cut.z = -2.0,
    min.z = -3.0,
    min.p = 1e-4,
    adj.p = 1.0
  )
  
  # Test p-values for scores in the ECDF range
  test_scores <- c(0, -0.5, -1.0, -1.5)
  pvals <- polylinkR:::.get_pvalue_lower(test_scores, 10, gpd_params, list(ecdf_fn))
  
  # P-values should be in [0, 1]
  expect_true(all(pvals >= 0 & pvals <= 1))
})

test_that(".estimate_setsize_covariance returns valid matrix", {
  # Test with a small set of sizes
  sizes <- c(5, 10, 20, 50)
  n.genes <- 1000
  
  K <- polylinkR:::.est_ss_cov(sizes, n.genes)
  
  # Should be a square matrix
  expect_equal(nrow(K), length(sizes))
  expect_equal(ncol(K), length(sizes))
  
  # Diagonal should be 1
  expect_equal(diag(K), rep(1, length(sizes)))
  
  # Should be symmetric
  expect_equal(K, t(K))
  
  # All values should be in [-1, 1]
  expect_true(all(K >= -1 & K <= 1))
})

test_that(".covariance_function returns expected values", {
  # Test exponential decay
  h <- c(0, 1, 2, 5, 10)
  ovlp <- c(1, 0.5, 0, 0, 0)
  
  result <- polylinkR:::.covariance_function(h, ovlp, beta1 = 1, beta2 = 0.5, scale = 2)
  
  # Should return numeric vector
  expect_type(result, "double")
  expect_length(result, length(h))
  
  # At h=0 with overlap=1, should return beta1 + beta2
  expect_equal(result[1], 1.5)
  
  # Should decay with distance
  expect_true(all(diff(result) < 0))
})

test_that(".cap_probabilities normalizes rows", {
  # Create a matrix with some values > maxP
  mm0 <- matrix(c(0.5, 0.3, 0.2,
                  0.6, 0.3, 0.1), nrow = 2, byrow = TRUE)
  
  result <- polylinkR:::.cap_probabilities(mm0, maxP = 0.4)
  
  # Each row should sum to 1 (within tolerance)
  expect_equal(rowSums(result), c(1, 1), tolerance = 1e-10)
  
  # No value should exceed maxP * rowSum
  expect_true(all(result <= 0.4))
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
