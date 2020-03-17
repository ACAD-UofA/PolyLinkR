test_that("Enrichment test and ploting works", {
  output = polylinkr(obj.info = Anatolia_EF_CLR, 
                     set.info = PolyLinkR_SetInfo,
                     set.obj = PolyLinkR_SetObj, 
                     n.cores = 4, 
                     emp.nruns = 10000, 
                     NN = 1000)
  
  expect_that(output, is_a("data.frame") )
  
})
