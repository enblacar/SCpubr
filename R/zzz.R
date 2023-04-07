# Add Start-Up message.
.onAttach <- function(...) {
  if (isFALSE(getOption("SCpubr.verbose"))){
    return()
  }
  
  # Print startup message.
  package_report(startup = TRUE)
}