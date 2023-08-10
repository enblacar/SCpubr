\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_LigandReceptorPlot", passive = TRUE)

  if (isTRUE(value)){
    liana_output <- readRDS(system.file("extdata/liana_output_example.rds", package = "SCpubr"))
    # Ligand Receptor analysis plot.
    p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
    p

  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

