# Add Start-Up message.
.onAttach <- function(...) {
      cli::cli_h1("SCpubr v{utils::packageVersion('SCpubr')}")
      
      cli::cli_alert_success(paste0(crayon_body("Have a look at extensive tutorials in "),
                                    crayon_key("SCpubr's book"),
                                    crayon_body(": https://enblacar.github.io/SCpubr-book/\n")))
      
      cli::cli_alert_success(paste0(crayon_body("If you use "), 
                                    crayon_key("SCpubr"),
                                    crayon_body(" in your research, please cite it accordingly: "),
                                    cli::style_italic("Blanco-Carmona, E. Generating publication ready visualizations for Single Cell transcriptomics using SCpubr. bioRxiv (2022)\ndoi:10.1101/2022.02.28.482303.")), wrap = TRUE)
      
      cli::cli_text("  ")
                             
      cli::cli_alert_info(paste0(crayon_body("If the package is useful to you, consider leaving a "),
                                   crayon_key("Star"),
                                   crayon_body(" in the GitHub repo: "),
                                   cli::style_italic("https://github.com/enblacar/SCpubr/stargazers\n")))
      
      cli::cli_alert_info(paste0(crayon_body("Keep track of the package "),
                                 crayon_key("updates"),
                                 crayon_body(" on Twitter ("),
                                 crayon_key("@Enblacar"),
                                 crayon_body(") or in: "),
                                 cli::style_italic("https://github.com/enblacar/SCpubr/blob/main/NEWS.md\n")))
      
      version <- utils::available.packages(repos = "http://cran.us.r-project.org")["SCpubr", "Version"]
       
      if (utils::packageVersion("SCpubr") < version){
         cli::cli_alert_warning(paste0(crayon_body("There is a new version of "),
                                       crayon_key("SCpubr (v"),
                                       crayon_key(version),
                                       crayon_key(")"),
                                       crayon_body(" available on "),
                                       crayon_key("CRAN"),
                                       crayon_body("! Please, consider updating the package to get the latest functionalities.\n")))
      }
      
      
      packages <- sort(unique(unlist(SCpubr::state_dependencies(return_dependencies = TRUE))))
      functions <- sort(unique(names(SCpubr::state_dependencies(return_dependencies = TRUE))))
      functions <- sapply(functions, SCpubr:::check_suggests, passive = TRUE)
      value <- c("Basics" = unname(functions[names(functions) == "core"]))
      functions <- c(value, functions[2:length(functions)])
      
      check_installed <- sapply(packages, requireNamespace, quietly = TRUE)
      max_length <- max(sapply(packages, nchar))
      max_length_functions <- max(sapply(names(functions), nchar))
      format_installed <- function(name, value, max_length){

        func_use <- ifelse(isTRUE(value), cli::col_green(cli::symbol$tick), cli::col_red(cli::symbol$cross))
        name_use <- ifelse(isTRUE(value), 
                           cli::ansi_align(cli::col_blue(name), max_length, align = "left"), 
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
        } else {
          print.vector <- c(print.vector, item)
          print.list[[item]] <- paste(print.vector, collapse = "     ")
          print.vector <- c()
        }
        
        if (counter == length(packages)){
          print.list[[item]] <- paste(print.vector, collapse = "     ")
          print.vector <- c()
        }
      }
      counter <- 0
      for(item in functions.use){
        counter <- counter + 1
        
        if (counter %% 4 != 0){
          print.vector.functions <- c(print.vector.functions, item)
        } else {
          print.vector.functions <- c(print.vector.functions, item)
          print.list.functions[[item]] <- paste(print.vector.functions, collapse = "     ")
          print.vector.functions <- c()
        }
        
        if (counter == length(functions)){
          print.list.functions[[item]] <- paste(print.vector.functions, collapse = "     ")
          print.vector.functions <- c()
        }
      }
      
      cli::cli_h2("Checking required packages")
      for (item in print.list){
        message(unname(item))
      }
      
      cli::cli_text("  ")
      
      cli::cli_h2("Checking available functions")
      for (item in print.list.functions){
        message(unname(item))
      }
      
      cli::cli_text("  ")
      
      cli::cli_alert_danger(crayon_body("To suppress this startup message, use: "))
      cli::cli_code('suppressMessages(library("SCpubr"))')
      message(cli::rule(col = "cyan"))
}