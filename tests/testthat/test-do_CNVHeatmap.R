if (base::isFALSE(dep_check[["do_CNVHeatmap"]])){

  testthat::test_that("do_BarPlot: CRAN essentials", {
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })

  testthat::test_that("do_BarPlot: PASS - normal cells all chromosomes", {

    testthat::skip_on_cran()

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations,
                                            flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                 infercnv_object = infercnv_object,
                                 using_metacells = FALSE,
                                 chromosome_locations = human_chr_locations,
                                 flip = TRUE,
                                 min.cutoff = 0.99,
                                 max.cutoff = 1.01)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                 infercnv_object = infercnv_object,
                                 using_metacells = FALSE,
                                 chromosome_locations = human_chr_locations,
                                 flip = TRUE,
                                 include_chr_arms = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                 infercnv_object = infercnv_object,
                                 using_metacells = FALSE,
                                 chromosome_locations = human_chr_locations,
                                 flip = TRUE,
                                 enforce_symmetry = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                 infercnv_object = infercnv_object,
                                 using_metacells = FALSE,
                                 chromosome_locations = human_chr_locations,
                                 flip = TRUE,
                                 values.show = TRUE,
                                 values.threshold = 1)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                 infercnv_object = infercnv_object,
                                 using_metacells = FALSE,
                                 chromosome_locations = human_chr_locations,
                                 flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations,
                                            flip = TRUE,
                                            group.by = c("seurat_clusters", "orig.ident", "annotation"))
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations,
                                            flip = FALSE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations,
                                            flip = FALSE,
                                            group.by = c("seurat_clusters", "orig.ident", "annotation"))
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations,
                                            flip = FALSE,
                                            return_object = TRUE)
    testthat::expect_type(out, "list")
    
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            group.by = c("seurat_clusters", "orig.ident"),
                                            chromosome_locations = human_chr_locations,
                                            flip = TRUE)
    testthat::expect_true(ggplot2::is_ggplot(out))
    
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            group.by = c("seurat_clusters", "orig.ident"),
                                            chromosome_locations = human_chr_locations,
                                            flip = FALSE,
                                            return_object = TRUE)
    testthat::expect_type(out, "list")
  })

  testthat::test_that("do_CNVHeatmap: PASS - normal cells one chromosome", {

    testthat::skip_on_cran()
    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object,
                                            using_metacells = FALSE,
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })

  testthat::test_that("do_CNVHeatmap: PASS - metacells all chromosomes", {
    testthat::skip_on_cran()

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })

  

  testthat::test_that("do_CNVHeatmap: PASS - group.by", {
    testthat::skip_on_cran()

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            group.by = "orig.ident",
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })


  testthat::test_that("do_CNVHeatmap: PASS - legend.position", {
    testthat::skip_on_cran()

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            legend.position = "right",
                                            legend.title = "test",
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            legend.position = "bottom",
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })

  testthat::test_that("do_CNVHeatmap: PASS - legend.position", {
    testthat::skip_on_cran()

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            legend.type = "normal",
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))

    out <- SCpubr::do_CNVHeatmap(sample = sample,
                                            infercnv_object = infercnv_object_metacells,
                                            using_metacells = TRUE,
                                            metacell_mapping = metacell_mapping,
                                            legend.type = "colorbar",
                                            chromosome_locations = human_chr_locations)
    testthat::expect_true(ggplot2::is_ggplot(out))
  })
}

