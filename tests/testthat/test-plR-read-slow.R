# Integration test: keep local and CRAN fast; skip on CI agents.
test_that("read_polylinkr_data succeeds on a small HumanPops extract", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()

  oi_full <- system.file("testthat", "humanpops_extract", "HumanPops_ObjInfo.txt", package = "polylinkR")
  si_full <- system.file("testthat", "humanpops_extract", "HumanPops_SetInfo.txt", package = "polylinkR")
  so_full <- system.file("testthat", "humanpops_extract", "HumanPops_SetObj.txt", package = "polylinkR")
  skip_if(oi_full == "", "Bundled HumanPops extract missing")

  d <- withr::local_tempdir()
  oi <- utils::read.delim(oi_full, nrows = 2000L, sep = "\t", stringsAsFactors = FALSE)
  utils::write.table(oi, file.path(d, "ObjInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  so <- utils::read.delim(so_full, sep = "\t", stringsAsFactors = FALSE)
  so <- so[so$objID %in% oi$objID, , drop = FALSE]
  # Find a set with >= 2 genes (the min(setID) set often has only 1 gene)
  set_counts <- table(so$setID)
  valid_sets <- names(set_counts[set_counts >= 2])
  stopifnot("No sets with >= 2 genes found" = length(valid_sets) > 0)
  so <- so[so$setID == valid_sets[1], , drop = FALSE]
  utils::write.table(so, file.path(d, "SetObj.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  si <- utils::read.delim(si_full, sep = "\t", stringsAsFactors = FALSE)
  si <- si[si$setID %in% unique(so$setID), , drop = FALSE]
  utils::write.table(si, file.path(d, "SetInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  out <- polylinkR::read_polylinkr_data(
    input_path = d,
    verbose = FALSE,
    min_set_size = 2L,
    max_set_size = Inf,
    merge_threshold = 0.95
  )
  expect_s3_class(out, "plR")
  expect_true(all(c("set.info", "obj.info", "set.obj") %in% names(out)))
})

# Backward compatibility test: ensure old function name still works
test_that("plR_read backward compatibility (deprecated) works", {
  testthat::skip_on_ci()
  testthat::skip_on_cran()
  
  # Create minimal test data
  d <- withr::local_tempdir()
  
  # Write minimal obj.info
  oi <- data.frame(
    objID = c("gene1", "gene2"),
    objStat = c(1.5, 2.3),
    chr = c("1", "1"),
    startpos = c(1000L, 2000L),
    endpos = c(1500L, 2500L),
    stringsAsFactors = FALSE
  )
  utils::write.table(oi, file.path(d, "ObjInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Write minimal set.info
  si <- data.frame(
    setID = c("set1"),
    stringsAsFactors = FALSE
  )
  utils::write.table(si, file.path(d, "SetInfo.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Write minimal set.obj
  so <- data.frame(
    setID = c("set1", "set1"),
    objID = c("gene1", "gene2"),
    stringsAsFactors = FALSE
  )
  utils::write.table(so, file.path(d, "SetObj.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Test that deprecated function still works but emits warning
  expect_warning(
    out <- polylinkR::plR_read(
      input.path = d,
      verbose = FALSE
    ),
    "deprecated"
  )
  expect_s3_class(out, "plR")
})