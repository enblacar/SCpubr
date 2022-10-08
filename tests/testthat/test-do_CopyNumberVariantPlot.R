
testthat::test_that("do_BarPlot: PASS - normal cells all chromosomes", {


  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE,
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - normal cells one chromosome", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object,
                                          using_metacells = FALSE,
                                          chromosome_focus = "2",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - metacells all chromosomes", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - metacells one chromosome", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - group.by", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          group.by = "orig.ident",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})


testthat::test_that("do_BarPlot: PASS - legend.position", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.position = "right",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.position = "bottom",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})

testthat::test_that("do_BarPlot: PASS - legend.position", {

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "normal",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "colorbar",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")

  out <- SCpubr::do_CopyNumberVariantPlot(sample = sample,
                                          infercnv_object = infercnv_object_metacells,
                                          using_metacells = TRUE,
                                          metacell_mapping = metacell_mapping,
                                          chromosome_focus = "2",
                                          legend.type = "colorsteps",
                                          chromosome_locations = human_chr_locations)
  testthat::expect_type(out, "list")
})
