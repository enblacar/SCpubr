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
    
    testthat::expect_warning({p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        assay = "SCT",
                                        slot = "data",
                                        nbin = 1,
                                        ctrl = 5)})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)
    testthat::expect_type(p, "list")
    
    testthat::expect_warning({p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "UCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"))
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5,
                                        use_viridis = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5,
                                        enforce_symmetry = TRUE)
    testthat::expect_type(p, "list")
    
  })
}


