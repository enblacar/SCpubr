# Add Start-Up message.
.onAttach <- function(...) {
  packageStartupMessage(crayon_body(paste(rep("-", 153), collapse = "")))
  packageStartupMessage(crayon::underline(crayon_key("SCpubr")))
  packageStartupMessage(paste0(crayon_body("\nIf you use "), 
                               crayon_key("SCpubr"),
                               crayon_body(" in your research, please cite it accordingly: \nBlanco-Carmona, E. Generating publication ready visualizations for Single Cell transcriptomics using SCpubr. bioRxiv (2022) doi:10.1101/2022.02.28.482303.\n")))
  packageStartupMessage(paste0(crayon_body("If the package is useful to you, consider leaving a "),
                               crayon_key("Star"),
                               crayon_body(" in the GitHub repo: https://github.com/enblacar/SCpubr/stargazers \n")))
  packageStartupMessage(paste0(crayon_body("Keep track of the package "),
                               crayon_key("updates"),
                               crayon_body(" on Twitter (@Enblacar) or in https://github.com/enblacar/SCpubr/blob/main/NEWS.md \n")))
  packageStartupMessage(crayon_body("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))"))
  packageStartupMessage(crayon_body(paste(rep("-", 153), collapse = "")))
}