\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_MetadataPlot", passive = TRUE)
  
  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/
    
    # Can also use a Seurat object.
    df <- data.frame(row.names = letters[1:5],
                     "A" = as.character(seq(1, 5)),
                     "B" = rev(as.character(seq(1, 5))))
    
    p <- SCpubr::do_MetadataPlot(from_df = TRUE,
                                 df = df)
    
  } else if (isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
