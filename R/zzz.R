# Add Start-Up message.
.onAttach <- function(...) {
  packageStartupMessage(paste(rep("-", 63), collapse = ""))
  packageStartupMessage("SCpubr")
  packageStartupMessage("\nIf you use SCpubr in your research, please cite it accordingly: \nBlanco-Carmona, E. Generating publication ready visualizations for Single Cell transcriptomics using SCpubr. bioRxiv (2022) doi:10.1101/2022.02.28.482303.\n")
  packageStartupMessage("Keep track of the package updates on Twitter (@Enblacar) or in https://github.com/enblacar/SCpubr/blob/main/NEWS.md \n")
  packageStartupMessage("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))")
  packageStartupMessage(paste(rep("-", 63), collapse = ""))
}
