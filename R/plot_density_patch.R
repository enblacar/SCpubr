#' Compatibility wrapper for Nebulosa with SeuratObject 5.0.0+
#'
#' Nebulosa is not compatible with SeuratObject >= 5.0.0 due to deprecated parameters.
#' This wrapper gracefully falls back to do_FeaturePlot() when needed.
#' 
#' @keywords internal
#' @noRd
.nebulosa_compat_wrapper <- function(object, features, reduction, dims, joint, slot = "data", 
                                     verbose = TRUE, ...) {
  
  # Check if we need the SeuratObject 5.0.0+ workaround
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    
    # Check for global option to suppress message (useful for testing)
    quiet_mode <- isTRUE(getOption("SCpubr.nebulosa.quiet", default = FALSE))
    
    if (isTRUE(verbose) && !quiet_mode) {
      message(paste0(add_info(), crayon_body("Nebulosa is not yet compatible with "),
                     crayon_key("SeuratObject >= 5.0.0"),
                     crayon_body(". Falling back to "),
                     crayon_key("do_FeaturePlot()"),
                     crayon_body(" instead.")))
    }
    
    # Return NULL to signal fallback is needed
    return(NULL)
    
  } else {
    # SeuratObject < 5.0.0, suppress warnings for aes_string()
    return(suppressWarnings({
      Nebulosa::plot_density(object = object,
                            features = features,
                            joint = joint,
                            reduction = reduction,
                            dims = dims)
    }))
  }
}
