# Add Start-Up message.
.onAttach <- function(...) {
  # nocov start
  if (base::isFALSE(getOption("SCpubr.verbose"))){
    return()
  }
  # nocov end
  
  # Print startup message.
  package_report(startup = TRUE,
                 extended = FALSE)
}