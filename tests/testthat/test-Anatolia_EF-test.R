test_that("Enrichment test", {
  out <- polylinkr(obj.info = Anatolia_EF_CLR, 
                     set.info = PolyLinkR_SetInfo,
                     set.obj = PolyLinkR_SetObj, 
                     n.cores = 4, 
                     emp.nruns = 100, 
                     NN = 1000)

   expect_output(str(out), "Classes ‘data.table’ and 'data.frame'")
   expect_equal(dim(out)[2], 6)
})
