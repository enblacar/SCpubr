sample <- SCpubr:::use_dataset()
metacell_mapping <- SCpubr:::metacell_mapping
infercnv_object <- SCpubr:::infercnv_object
infercnv_object_metacells <- SCpubr:::infercnv_object_metacells


testthat::test_that("do_BarPlot: PASS - normal cells all chromosomes", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - normal cells one chromosome", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE,
                                          chromosome_focus = "2")
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - metacells all chromosomes", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - metacells one chromosome", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2")
  testthat::expect_type(out, "list")
})


testthat::test_that("do_BarPlot: PASS - legend.position", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.position = "right")
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.position = "bottom")
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - legend.position", {
  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "normal")
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "colorbar")
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "colorsteps")
  testthat::expect_type(out, "list")
})
