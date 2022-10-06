sample <- SCpubr:::use_dataset()
testthat::test_that("do_RidgePlot: PASS - default", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_RidgePlot: PASS - continuous scale", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            continuous_scale = TRUE,
                            viridis_direction = 1)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            continuous_scale = TRUE,
                            viridis_direction = -1)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_RidgePlot: PASS - group.by", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            group.by = "orig.ident")
  testthat::expect_type(p, "list")
})

testthat::test_that("do_RidgePlot: PASS - quantiles", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            continuous_scale = TRUE,
                            compute_custom_quantiles = TRUE,
                            compute_quantiles = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            continuous_scale = TRUE,
                            compute_custom_quantiles = TRUE,
                            compute_quantiles = TRUE,
                            quantiles = c(0.1, 0.5, 0.9))
  testthat::expect_type(p, "list")
})

testthat::test_that("do_RidgePlot: PASS - distribution tails", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nFeature_RNA",
                            continuous_scale = TRUE,
                            compute_quantiles = TRUE,
                            compute_distribution_tails = TRUE)
  testthat::expect_type(p, "list")

  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nCount_RNA",
                            continuous_scale = TRUE,
                            compute_quantiles = TRUE,
                            compute_distribution_tails = TRUE,
                            prob_tails = 0.4)
  testthat::expect_type(p, "list")
})

testthat::test_that("do_RidgePlot: PASS - distribution tails", {
  p <- SCpubr::do_RidgePlot(sample = sample,
                            feature = "nFeature_RNA",
                            continuous_scale = TRUE,
                            compute_quantiles = TRUE,
                            color_by_probabilities = TRUE)
  testthat::expect_type(p, "list")
})


