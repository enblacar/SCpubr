sample <- use_dataset()
p <- do_DimPlot(sample)
figure_path <- "./"
testthat::test_that("save_Plot: PASS - all", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "all"))
})

testthat::test_that("save_Plot: PASS - publication", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "publication"))
})

testthat::test_that("save_Plot: PASS - jpeg", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "jpeg"))
})

testthat::test_that("save_Plot: PASS - png", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "png"))
})

testthat::test_that("save_Plot: PASS - pdf", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "pdf"))
})

testthat::test_that("save_Plot: PASS - tiff", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "tiff"))
})

testthat::test_that("save_Plot: PASS - svg", {
  testthat::expect_silent(save_Plot(plot = p,
                                    figure_path = figure_path,
                                    file_name = "test",
                                    output_format = "svg"))
})

unlink(paste0(figure_path, "test.svg"))
unlink(paste0(figure_path, "test.jpeg"))
unlink(paste0(figure_path, "test.pdf"))
unlink(paste0(figure_path, "test.tiff"))
unlink(paste0(figure_path, "test.png"))
