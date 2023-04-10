if (isFALSE(dep_check[["do_DiffusionMapPlot"]])){
  
  testthat::test_that("do_DiffusionMapPlot: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_DiffusionMapPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     verbose = FALSE)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_DiffusionMapPlot: PASS - default", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_DiffusionMapPlot(sample = sample,
                                     input_gene_list =  genes,
                                     subsample = NA,
                                     nbin = 1,
                                     ctrl = 5,
                                     reduction = "umap",
                                     dims = 1:2,
                                     return_object = TRUE,
                                     verbose = FALSE)
    testthat::expect_type(p, "list")
  })
}


