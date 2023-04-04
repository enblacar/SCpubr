# Add Start-Up message.
.onAttach <- function(...) {
  if (isFALSE(getOption("SCpubr.verbose"))){
    return()
  }
  
  # Check if cli is installed, otherwise just output normal text and return..
  if (isFALSE(requireNamespace("cli", quietly = TRUE))){
    packageStartupMessage(paste(rep("-", 63), collapse = ""))
    packageStartupMessage('This is a placeholder message. Please install "cli" package to have an optimal experience using the package.')
    packageStartupMessage(paste(rep("-", 63), collapse = ""))
    packageStartupMessage("\n\n\nSCpubr")
    packageStartupMessage("\nIf you use SCpubr in your research, please cite it accordingly: \nBlanco-Carmona, E. Generating publication ready visualizations for Single Cell transcriptomics using SCpubr. bioRxiv (2022) doi:10.1101/2022.02.28.482303.\n")
    packageStartupMessage("If the package is useful to you, consider leaving a Star in the GitHub repo: https://github.com/enblacar/SCpubr/stargazers \n")
    packageStartupMessage("Keep track of the package updates on Twitter (@Enblacar) or in https://github.com/enblacar/SCpubr/blob/main/NEWS.md \n")
    packageStartupMessage("To suppress this startup message, use: \nsuppressPackageStartupMessages(library('SCpubr'))")
    packageStartupMessage(paste(rep("-", 63), collapse = ""))
    
    return()
  }
  tip_rule <- cli::rule(left = "General", width = nchar("General") + 6)
  
  tutorials <- paste0(add_info(),
                      crayon_body("Have a look at extensive tutorials in "),
                      crayon_key(cli::style_hyperlink(text = "SCpubr's book",
                                                      url = "https://enblacar.github.io/SCpubr-book/")),
                      crayon_body("."))
  
  cite <- paste0(add_tick(),
                 crayon_body("If you use "),
                 crayon_key("SCpubr"),
                 crayon_body(" in your research, please "),
                 crayon_key(cli::style_hyperlink(text = "cite it accordingly",
                                                 url = "https://www.biorxiv.org/content/10.1101/2022.02.28.482303v1")),
                 crayon_body("."))
  
  stars <- paste0(add_star(),
                  crayon_body("If the package is useful to you, consider leaving a "),
                  crayon_key("Star"),
                  crayon_body(" in the "),
                  crayon_key(cli::style_hyperlink(text = "GitHub repository",
                                                  url = "https://github.com/enblacar/SCpubr")),
                  crayon_body("."))
  
  updates <- paste0(cli::style_bold(cli::col_blue("!")),
                    crayon_body(" Keep track of the package "),
                    crayon_key("updates"),
                    crayon_body(" on Twitter ("),
                    crayon_key(cli::style_hyperlink(text = "@Enblacar",
                                                    url = "https://twitter.com/Enblacar")),
                    crayon_body(") or in the "),
                    crayon_key(cli::style_hyperlink(text = "Official NEWS website",
                                                    url = "https://github.com/enblacar/SCpubr/blob/main/NEWS.md")),
                    crayon_body("."))
  
  plotting <- paste0(cli::style_bold(cli::col_red(cli::symbol$heart)), " ", crayon_body("Happy plotting!"))
  
  updates_check <- cli::rule(left = "Checking package updates", width = nchar("Checking package updates") + 6)
  
  cran_version <- utils::available.packages(repos = "http://cran.us.r-project.org")["SCpubr", "Version"]
  system_version <- as.character(utils::packageVersion("SCpubr"))
  
  l1 <- nchar(system_version)
  l2 <- nchar(cran_version)
  max_length <- max(c(l1, l2))
  
  if (rev(strsplit(as.character(system_version), split = "\\.")[[1]])[1] >= 9000){
    parts <- strsplit(as.character(system_version), split = "\\.")[[1]]
    parts[[length(parts)]] <- cli::col_yellow(parts[[length(parts)]])
    system_version <- paste(parts, collapse = ".")
  }
  
  header <- cli::rule(left = paste0(crayon_body("SCpubr "),
                                    crayon_key(system_version)), line_col = "cadetblue")
  
  system_version_message <- paste0(cli::col_magenta("System: "), cli::ansi_align(crayon_key(system_version), max_length, align = "right"))
  cran_version_message <- paste0(cli::col_magenta("CRAN:   "), cli::ansi_align(crayon_key(cran_version), max_length, align = "right"))
  
  
  
  if (system_version < cran_version){
    veredict_message <- paste0(cli::col_yellow(cli::symbol$warning,
                                               crayon_body(" There is a "),
                                               crayon_key("new version"),
                                               crayon_body(" available on"),
                                               crayon_key("CRAN"),
                                               crayon_body("!")))
  } else {
    veredict_message <- paste0(cli::col_green(cli::symbol$tick,
                                              crayon_body(" Installation is "),
                                              crayon_key("up to date"),
                                              crayon_body(" with latest "),
                                              crayon_key("CRAN"),
                                              crayon_body(" version!")))
  }
  
  
  
  packages <- sort(unique(unlist(check_dependencies(return_dependencies = TRUE))))
  functions <- sort(unique(names(check_dependencies(return_dependencies = TRUE))))
  
  if (rev(strsplit(as.character( as.character(utils::packageVersion("SCpubr"))), split = "\\.")[[1]])[1] >= 9000){
    names.use <- unname(sapply(functions, function(x){if (x %in% c("do_SankeyPlot", "do_PseudotimePlot", "do_LigandReceptorPlot", "save_Plot")){x <- paste0(x, cli::col_yellow(" | DEV"))} else {x}}))
    functions <- sapply(functions, check_suggests, passive = TRUE)
    names(functions) <- names.use
  } else {
    functions <- functions[!(functions %in% c("do_SankeyPlot", "do_PseudotimePlot", "do_LigandReceptorPlot", "save_Plot"))]
    functions <- sapply(functions, check_suggests, passive = TRUE)
  }
  
  
  functions <- functions[names(functions) != "Essentials"]
  
  check_installed <- sapply(packages, requireNamespace, quietly = TRUE)
  max_length <- max(sapply(packages, nchar))
  max_length_functions <- max(sapply(names(functions), nchar))
  format_installed <- function(name, value, max_length){
    func_use <- ifelse(isTRUE(value), cli::col_green(cli::symbol$tick), cli::col_red(cli::symbol$cross))
    name_use <- ifelse(isTRUE(value),
                       cli::ansi_align(crayon_key(name), max_length, align = "left"),
                       cli::ansi_align(cli::col_red(name), max_length, align = "left"))
    paste0(func_use, " ", name_use)
  }
  
  packages <- c()
  for(item in names(check_installed)){
    packages <- c(packages, format_installed(name = item, value = check_installed[[item]], max_length = max_length))
  }
  
  functions.use <- c()
  for(item in names(functions)){
    functions.use <- c(functions.use, format_installed(name = item, value = functions[[item]], max_length = max_length_functions))
  }
  
  counter <- 0
  print.list <- list()
  print.list.functions <- list()
  print.vector <- c()
  print.vector.functions <- c()
  for(item in packages){
    counter <- counter + 1
    
    if (counter %% 5 != 0){
      print.vector <- c(print.vector, item)
      if (counter == length(packages)){
        print.list[[item]] <- paste(print.vector, collapse = "     ")
        print.vector <- c()
      }
    } else {
      print.vector <- c(print.vector, item)
      print.list[[item]] <- paste(print.vector, collapse = "     ")
      print.vector <- c()
    }
  }
  
  counter <- 0
  for(item in functions.use){
    counter <- counter + 1
    
    if (counter %% 4 != 0){
      print.vector.functions <- c(print.vector.functions, item)
      if (counter == length(functions.use)){
        print.list.functions[[item]] <- paste(print.vector.functions, collapse = "     ")
      }
    } else {
      print.vector.functions <- c(print.vector.functions, item)
      print.list.functions[[item]] <- paste(print.vector.functions, collapse = "     ")
      print.vector.functions <- c()
    }
    
    
  }
  
  packages_check <- cli::rule(left = "Required packages", width = nchar("Required packages") + 6)
  functions_check <- cli::rule(left = "Available functions", width = nchar("Available functions") + 6)
  
  functions_tip <- paste0(cli::style_bold(cli::col_cyan(cli::symbol$info)), 
                          crayon_body(" Check the package requirements function-wise with: "), 
                          cli::style_italic(crayon_key('SCpubr::check_dependencies()')))
  
  tip_rule <- cli::rule(left = "Tips!", width = nchar("Tips!") + 6)
  
  tip_message <- paste0(cli::style_bold(cli::col_cyan(cli::symbol$info)), 
                        crayon_body(" To adjust package messages to dark mode themes, use: "), 
                        cli::style_italic(crayon_key('options("SCpubr.darkmode" = TRUE)')))
  
  disable_message <- paste0(cli::style_bold(cli::col_red(cli::symbol$cross)), 
                            crayon_body(" To suppress this startup message, use: "), 
                            cli::style_italic(crayon_key('suppressPackageStartupMessages(library(SCpubr))\n')),
                            cli::style_bold(cli::col_red(cli::symbol$cross)), 
                            crayon_body(" Alternatively, you can also set the following option: "),
                            cli::style_italic(crayon_key('options("SCpubr.verbose" = FALSE)\n')),
                            crayon_body("  And then load the package normally (and faster) as: "),
                            cli::style_italic(crayon_key('library(SCpubr)')))
  
  end_rule <- cli::rule(col = "cadetblue")
  
  # Mount all individual messages into a big one that will be then be printed as a packageStartupMessage.
  msg_wrap <- paste0("\n", "\n", 
                     header, "\n", "\n",
                     tutorials, "\n", "\n",
                     cite, "\n", "\n",
                     stars, "\n", "\n",
                     updates, "\n", "\n",
                     plotting, "\n", "\n", "\n",
                     updates_check, "\n", "\n",
                     cran_version_message, "\n",
                     system_version_message, "\n", "\n",
                     veredict_message, "\n", "\n", "\n",
                     packages_check, "\n", "\n",
                     paste(print.list, collapse = "\n"), "\n", "\n", "\n",
                     functions_check, "\n", "\n",
                     paste(print.list.functions, collapse = "\n"), "\n","\n",
                     functions_tip, "\n", "\n", "\n",
                     tip_rule, "\n", "\n",
                     tip_message, "\n", "\n",
                     disable_message, "\n", "\n",
                     end_rule)
  
  # This allows for the package to be silenced on demand.
  rlang::inform(msg_wrap, class = "packageStartupMessage")
}