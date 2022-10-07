

testthat::test_that("save_Plot: PASS - no file", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            output_format = "svg"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - no file path", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            file_name = "test",
                                            output_format = "svg"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - null file path", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            file_name = "test",
                                            output_format = "svg"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - no file path", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = paste0(figure_path, "/deleteme"),
                                            file_name = "test",
                                            output_format = "svg"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: FAIL - wrong output format", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_error(SCpubr::save_Plot(plot = p,
                                           figure_path = figure_path,
                                           file_name = "test",
                                           output_format = "wrong"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: FAIL - all and publication at the same time.", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_error(SCpubr::save_Plot(plot = p,
                                           figure_path = figure_path,
                                           file_name = "test",
                                           output_format = c("all", "publication")))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - all", {

  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "all"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "all"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "all"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.chord,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "all"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - publication", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "publication"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "publication"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "publication"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.chord,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "publication"))

  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)

})

testthat::test_that("save_Plot: PASS - jpeg", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "jpeg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "jpeg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "jpeg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.chord,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "jpeg"))
  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - png", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "png"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "png"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "png"))
})

testthat::test_that("save_Plot: PASS - pdf", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "pdf"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "pdf"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "pdf"))
  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - tiff", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "tiff"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "tiff"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "tiff"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.chord,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "tiff"))
  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})

testthat::test_that("save_Plot: PASS - svg", {
  sample <- SCpubr:::test_list$sample
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()

  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "svg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.heatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "svg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.pheatmap,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "svg"))

  testthat::expect_silent(SCpubr::save_Plot(plot = p.chord,
                                            figure_path = figure_path,
                                            file_name = "test",
                                            output_format = "svg"))
  unlink(paste0(figure_path, "*.svg"))
  unlink(paste0(figure_path, "test.jpeg"))
  unlink(paste0(figure_path, "test.pdf"))
  unlink(paste0(figure_path, "test.tiff"))
  unlink(paste0(figure_path, "test.png"))
  unlink(paste0(figure_path, "/deleteme"), recursive = T)
})




