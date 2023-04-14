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
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:11],
                  "C" = rownames(sample)[12:19])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE,
                                         compute_robustness = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = TRUE,
                                         compute_robustness = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = FALSE,
                                         compute_robustness = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = TRUE,
                                         compute_robustness = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = FALSE,
                                         compute_robustness = TRUE)
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE,
                                         return_object = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - robustness", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE,
                                         compute_robustness = TRUE)
    testthat::expect_type(p, "list")
    
    suppressMessages({testthat::expect_message({p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = TRUE,
                                                                   flip = TRUE,
                                                                   compute_robustness = TRUE)})})
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:3],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[9:15])
    testthat::expect_error({SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                                            input_gene_list =  genes,
                                                            subsample = 100,
                                                            nbin = 1,
                                                            ctrl = 5,
                                                            verbose = FALSE,
                                                            compute_robustness = TRUE)})
    genes <- list("A" = rownames(sample)[1:15],
                  "B" = rownames(sample)[16:40],
                  "C" = rownames(sample)[41:80])
    
    SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE,
                                    compute_robustness = FALSE)
    testthat::expect_type(p, "list")
    
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - symmetry", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = FALSE,
                                         use_viridis = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = FALSE,
                                         use_viridis = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = TRUE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - add enrichment", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         add.enrichment = TRUE,
                                         use_viridis = TRUE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         add.enrichment = TRUE,
                                         use_viridis = TRUE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    suppressMessages({testthat::expect_message({ p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                                                    input_gene_list =  genes,
                                                                    subsample = 100,
                                                                    nbin = 1,
                                                                    ctrl = 5,
                                                                    verbose = TRUE,
                                                                    add.enrichment = TRUE,
                                                                    use_viridis = TRUE)})})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         add.enrichment = TRUE,
                                         use_viridis = TRUE,
                                         flavor = "UCell")
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         add.enrichment = TRUE,
                                         use_viridis = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = TRUE,
                                         add.enrichment = TRUE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - flip", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - cutoffs", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         min.cutoff = -0.25,
                                         max.cutoff = 0.25)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - multiple group.by", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - verbose", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    testthat::expect_message({p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = TRUE)})
    
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - underscores", {
    testthat::skip_on_cran()
    genes <- list("_A" = rownames(sample)[1:5],
                  "_B" = rownames(sample)[6:10],
                  "_C" = rownames(sample)[11:15])
    
    testthat::expect_warning({p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = FALSE)}) 
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_AffinityAnalysisPlot: PASS - different length of gene sets", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:15],
                  "C" = rownames(sample)[15:30])
    
    p <- SCpubr::do_AffinityAnalysisPlot(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
  })
  
}


