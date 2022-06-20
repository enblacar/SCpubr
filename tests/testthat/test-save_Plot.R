sample <- SCpubr:::use_dataset()
p <- SCpubr::do_DimPlot(sample)
p.heatmap <- SCpubr::do_CorrelationPlot(sample)
data <- p.heatmap@ht_list$`Pearson coef.`@matrix
p.pheatmap <- pheatmap::pheatmap(data)
figure_path <- "./"

testthat::test_that("save_Plot: PASS - no file", {
  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = figure_path,
                                            output_format = "svg"))
})

testthat::test_that("save_Plot: PASS - no file path", {
  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            file_name = "test",
                                            output_format = "svg"))
})

testthat::test_that("save_Plot: PASS - no file path", {
  testthat::expect_silent(SCpubr::save_Plot(plot = p,
                                            figure_path = paste0(figure_path, "/deleteme"),
                                            file_name = "test",
                                            output_format = "svg"))
})

testthat::test_that("save_Plot: FAIL - wrong output format", {
  testthat::expect_error(SCpubr::save_Plot(plot = p,
                                           figure_path = figure_path,
                                           file_name = "test",
                                           output_format = "wrong"))
})

testthat::test_that("save_Plot: PASS - all", {
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
})

testthat::test_that("save_Plot: PASS - publication", {
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

})

testthat::test_that("save_Plot: PASS - jpeg", {
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
})

testthat::test_that("save_Plot: PASS - png", {
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
})

testthat::test_that("save_Plot: PASS - tiff", {
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
})

testthat::test_that("save_Plot: PASS - svg", {
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
})



unlink(paste0(figure_path, "*.svg"))
unlink(paste0(figure_path, "test.jpeg"))
unlink(paste0(figure_path, "test.pdf"))
unlink(paste0(figure_path, "test.tiff"))
unlink(paste0(figure_path, "test.png"))
unlink(paste0(figure_path, "/deleteme"), recursive = T)
