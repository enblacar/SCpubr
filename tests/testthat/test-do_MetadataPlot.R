if (isFALSE(dep_check[["do_MetadataPlot"]])){
  
  testthat::test_that("do_MetadataPlot: CRAN essentials", {
    df <- data.frame(row.names = letters[1:5],
                     "A" = as.character(seq(1, 5)),
                     "B" = rev(as.character(seq(1, 5))))
    
    p <- SCpubr::do_MetadataPlot(from_df = TRUE,
                                 df = df)
    testthat::expect_type(p, "list")
    
    
  })
  
  testthat::test_that("do_MetadataPlot: PASS - default", {
    testthat::skip_on_cran()
    
    df <- data.frame(row.names = letters[1:5],
                     "A" = as.character(seq(1, 5)),
                     "B" = rev(as.character(seq(1, 5))))
    
    p <- SCpubr::do_MetadataPlot(from_df = TRUE,
                                 df = df)
    testthat::expect_type(p, "list")
  })
}