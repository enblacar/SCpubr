if (isFALSE(dep_check[["do_AffinityAnalysisPlot"]])){
  
  testthat::test_that("do_AffinityAnalysisPlot: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = NA,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - default", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = NA,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
  })
}


