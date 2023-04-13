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
                     "B" = rev(as.character(seq(1, 5))),
                     "C" = c("1", "2", "3", "5", "7"))
    
    p <- SCpubr::do_MetadataPlot(from_df = TRUE,
                                 df = df,
                                 flip = FALSE,
                                 legend.symbol.size = 2)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_MetadataPlot(from_df = TRUE,
                                 df = df,
                                 flip = TRUE)
    testthat::expect_type(p, "list")
    
    sample$labelling <- sample(c("A", "B"), ncol(sample), replace = TRUE)
    p <- SCpubr::do_MetadataPlot(sample = sample,
                                 group.by = "labelling",
                                 metadata = "orig.ident",
                                 flip = FALSE)
    testthat::expect_type(p, "list")
    
    p <- SCpubr::do_MetadataPlot(sample = sample,
                                 group.by = "labelling",
                                 metadata = "orig.ident",
                                 flip = TRUE)
    testthat::expect_type(p, "list")
  })
}