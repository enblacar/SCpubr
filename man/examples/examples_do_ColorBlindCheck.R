\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_ColorBlindCheck", passive = TRUE)
  
  if (isTRUE(value)){
    # Generate a color wheel based on a single value.
    colors <- c("red", "green", "blue")
    p <- SCpubr::do_ColorBlindCheck(colors.use = colors)
    
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}
