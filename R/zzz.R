# Add Start-Up message.
.onAttach <- function(...) {
  if(!requireNamespace("crayon", quietly = T)){
    packageStartupMessage(paste(rep("-", 63), collapse = ""))
    packageStartupMessage("SCpubr")
    packageStartupMessage("\nIf you use SCpubr in your research, please cite it accordingly: \nutils::citation(package = 'SCpubr')\n")
    packageStartupMessage("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))")
    packageStartupMessage(paste(rep("-", 63), collapse = ""))
  } else {
    packageStartupMessage(crayon::yellow(crayon::bold(paste(rep("-", 63), collapse = ""))))
    packageStartupMessage(crayon::yellow("SCpubr"))
    packageStartupMessage(crayon::yellow("\nIf you use SCpubr in your research, please cite it accordingly: \nutils::citation(package = 'SCpubr')\n"))
    packageStartupMessage(crayon::green("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))"))
    packageStartupMessage(crayon::yellow(crayon::bold(paste(rep("-", 63), collapse = ""))))
  }
}
