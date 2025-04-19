if (base::isFALSE(dep_check[["do_ActivityHeatmap"]])){
  
  testthat::test_that("do_ActivityHeatmap: CRAN essentials", {
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = NA,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = NA,
                                    nbin = 1,
                                    ctrl = 5,
                                    enforce_symmetry = FALSE,
                                    verbose = FALSE)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - default", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:11],
                  "C" = rownames(sample)[12:24])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = TRUE,
                                    statistic = "wmean",
                                    flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    group.by = "orig.ident",
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE,
                                    flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    min.cutoff = 0,
                                    max.cutoff = 0.1,
                                    verbose = FALSE,
                                    flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE,
                                    values.show = TRUE,
                                    values.threshold = 0.1,
                                    enforce_symmetry = TRUE,
                                    flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE,
                                    values.show = TRUE,
                                    values.threshold = 0.1,
                                    enforce_symmetry = FALSE,
                                    flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE,
                                    flip = FALSE,
                                    values.show = TRUE,
                                    values.threshold = 0.2,
                                    return_object = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE,
                                         return_object = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - robustness", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    suppressMessages({testthat::expect_message({p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = TRUE,
                                                                   flip = TRUE)})})
    testthat::expect_type(p, "list")
    
    genes <- list("A" = rownames(sample)[1:3],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[9:15])
    testthat::expect_error({SCpubr::do_ActivityHeatmap(sample = sample,
                                                            input_gene_list =  genes,
                                                            subsample = 100,
                                                            nbin = 1,
                                                            ctrl = 5,
                                                            verbose = FALSE)})
    genes <- list("A" = rownames(sample)[1:15],
                  "B" = rownames(sample)[16:40],
                  "C" = rownames(sample)[41:80])
    
    SCpubr::do_ActivityHeatmap(sample = sample,
                                    input_gene_list =  genes,
                                    subsample = 100,
                                    nbin = 1,
                                    ctrl = 5,
                                    verbose = FALSE)
    testthat::expect_type(p, "list")
    
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - symmetry", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = FALSE,
                                         use_viridis = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = FALSE,
                                         use_viridis = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = TRUE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - add enrichment", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         use_viridis = TRUE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         use_viridis = TRUE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
    
    suppressMessages({testthat::expect_message({ p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                                                    input_gene_list =  genes,
                                                                    subsample = 100,
                                                                    nbin = 1,
                                                                    ctrl = 5,
                                                                    verbose = TRUE,
                                                                    use_viridis = TRUE)})})
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         use_viridis = TRUE,
                                         flavor = "UCell")
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         use_viridis = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         enforce_symmetry = TRUE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - flip", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         flip = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - cutoffs", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE,
                                         min.cutoff = -0.25,
                                         max.cutoff = 0.25)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - multiple group.by", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         group.by = c("seurat_clusters", "orig.ident"),
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - verbose", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:10],
                  "C" = rownames(sample)[11:15])
    
    testthat::expect_message({p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = TRUE)})
    
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - underscores", {
    testthat::skip_on_cran()
    genes <- list("_A" = rownames(sample)[1:5],
                  "_B" = rownames(sample)[6:10],
                  "_C" = rownames(sample)[11:15])
    
    testthat::expect_warning({p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                                                   input_gene_list =  genes,
                                                                   subsample = 100,
                                                                   nbin = 1,
                                                                   ctrl = 5,
                                                                   verbose = FALSE)}) 
    testthat::expect_type(p, "list")
  })
  
  testthat::test_that("do_ActivityHeatmap: PASS - different length of gene sets", {
    testthat::skip_on_cran()
    genes <- list("A" = rownames(sample)[1:5],
                  "B" = rownames(sample)[6:15],
                  "C" = rownames(sample)[15:30])
    
    p <- SCpubr::do_ActivityHeatmap(sample = sample,
                                         input_gene_list =  genes,
                                         subsample = 100,
                                         nbin = 1,
                                         ctrl = 5,
                                         verbose = FALSE)
    testthat::expect_type(p, "list")
  })
  
}


