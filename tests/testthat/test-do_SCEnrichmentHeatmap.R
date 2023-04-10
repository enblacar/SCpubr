if (isFALSE(dep_check[["do_SCEnrichmentHeatmap"]])){
  
  testthat::test_that("do_SCEnrichmentHeatmap: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_SCEnrichmentHeatmap: PASS - default", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)
    testthat::expect_type(p, "list")
  })
}


