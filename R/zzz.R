# Add Start-Up message.
.onAttach <- function(...) {
  packageStartupMessage(paste(rep("-", 63), collapse = ""))
  packageStartupMessage("SCpubr")
  packageStartupMessage("\nIf you use SCpubr in your research, please cite it accordingly: \nutils::citation(package = 'SCpubr')\n")
  packageStartupMessage("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))")
  packageStartupMessage(paste(rep("-", 63), collapse = ""))
}
