if(isFALSE(dep_check[["do_LigandReceptorPlot"]])){
  testthat::test_that("do_LigandReceptorPlot: CRAN essentials", {
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output", {
    testthat::skip_on_cran()
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = TRUE, viridis.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = TRUE, viridis.direction = -1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = FALSE, sequential.direction = -1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = TRUE, use_viridis = FALSE, sequential.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = TRUE, viridis.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = TRUE, viridis.direction = -1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = FALSE, sequential.direction = -1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, dot_border = FALSE, use_viridis = FALSE, sequential.direction = 1)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = TRUE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, plot.grid = TRUE, dot_border = FALSE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output, 
                                       keep_source = c("NK", "B"),
                                       keep_target = c("CD8 T"))
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output different n", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50)
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_LigandReceptorPlot: PASS - split.by", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       split.by = "ligand.complex")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       split.by = "receptor.complex")
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: PASS - from output, angle ", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 0)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 45)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       axis.text.x.angle = 90)
    testthat::expect_type(p, "list")
  })


  testthat::test_that("do_LigandReceptorPlot: PASS - from output legend.position", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       legend.position = "bottom")
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       legend.position = "right")
    testthat::expect_type(p, "list")
  })


  

  testthat::test_that("do_LigandReceptorPlot: PASS - sort interactions", {
    testthat::skip_on_cran()
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  TRUE)
    testthat::expect_type(p, "list")

    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 50,
                                       sort_interactions_alphabetically =  FALSE)
    testthat::expect_type(p, "list")
  })

  testthat::test_that("do_LigandReceptorPlot: FAIL - wrong parameters", {
    testthat::skip_on_cran()
    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          font.type = "wrong")})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          legend.type = "wrong")})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          axis.text.x.angle = 10)})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          font.type = "wrong")})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          legend.position = "wrong")})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          grid.type = "wrong")})

    testthat::expect_error({SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                                          split.by = "wrong")})

  })
}

