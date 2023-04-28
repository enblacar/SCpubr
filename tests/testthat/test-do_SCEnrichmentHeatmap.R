if (base::isFALSE(dep_check[["do_SCEnrichmentHeatmap"]])){
  
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
                                        ctrl = 5,
                                        cluster = FALSE,
                                        features.order = c("B", "C", "A"))
    testthat::expect_type(p, "list")
    
    genes <- list("A_A" = rownames(sample)[1:5],
                  "B_A" = rownames(sample)[6:10],
                  "C_A" = rownames(sample)[11:15])
    
    suppressWarnings({testthat::expect_warning({p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                                                  input_gene_list = genes,
                                                                  flavor = "Seurat",
                                                                  assay = "SCT",
                                                                  slot = "data",
                                                                  nbin = 1,
                                                                  ctrl = 5)})})
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    testthat::expect_error({p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                                                  input_gene_list = "EPC1",
                                                                  flavor = "Seurat",
                                                                  assay = "SCT",
                                                                  nbin = 1,
                                                                  ctrl = 5)})
    
    sample$test <- as.factor(sample$seurat_clusters)
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        group.by = "test",
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5])
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        subsample = 100,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5)
    testthat::expect_type(p, "list")
    
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5,
                                        cluster = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5,
                                        cluster = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("orig.ident", "seurat_clusters"),
                                        metadata.colors = list("orig.ident" = c("Cell" = "red")))
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5,
                                        proportional.size = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "AUCell",
                                        assay = "SCT",
                                        nbin = 1,
                                        ctrl = 5,
                                        proportional.size = FALSE)
    testthat::expect_type(p, "list")
   
    
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
                                        use_viridis = TRUE,
                                        viridis.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5,
                                        use_viridis = TRUE,
                                        viridis.direction = -1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5,
                                        use_viridis = FALSE,
                                        sequential.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        metadata = c("seurat_clusters", "orig.ident"),
                                        min.cutoff = 0,
                                        max.cutoff = 0.5,
                                        use_viridis = FALSE,
                                        sequential.direction = -1)
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
    
    p <- SCpubr::do_SCEnrichmentHeatmap(sample = sample,
                                        input_gene_list = genes,
                                        flavor = "Seurat",
                                        nbin = 1,
                                        ctrl = 5,
                                        return_object = TRUE)
    testthat::expect_type(p, "list")
    
  })
}


