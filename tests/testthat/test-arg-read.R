test_that(".arg_check read rejects set.merge outside (0, 1]", {
  e <- read_arg_env(set.merge = 0L)
  expect_error(polylinkR:::.arg_check("read", e), "set.merge")
  e2 <- read_arg_env(set.merge = 1.5)
  expect_error(polylinkR:::.arg_check("read", e2), "set.merge")
})

test_that(".arg_check read rejects min.set.n and max.set.n inconsistency", {
  e <- read_arg_env(min.set.n = 5L, max.set.n = 3L)
  expect_error(polylinkR:::.arg_check("read", e), "min.set.n > max.set.n")
})

test_that(".arg_check read rejects both set.in and set.out", {
  e <- read_arg_env(set.in = 1L, set.out = 2L)
  expect_error(polylinkR:::.arg_check("read", e), "set.in or set.out")
})

test_that(".arg_check read rejects invalid obj.stat.fun", {
  e <- read_arg_env(obj.stat.fun = "invalid")
  expect_error(polylinkR:::.arg_check("read", e), "obj.stat.fun")
})

test_that(".arg_check read rejects invalid map.fun", {
  e <- read_arg_env(map.fun = "mercator")
  expect_error(polylinkR:::.arg_check("read", e), "kosambi")
})

test_that(".arg_check read rejects bin.size out of range", {
  e <- read_arg_env(bin.size = 10L)
  expect_error(polylinkR:::.arg_check("read", e), "bin.size")
})
