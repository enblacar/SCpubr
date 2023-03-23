# Add Start-Up message.
.onAttach <- function(...) {
      cli::cli_h1("SCpubr v{utils::packageVersion('SCpubr')}")
      
      cli::cli_text(" ")
  
      cli::cli_alert_info(paste0(crayon_body("Have a look at extensive tutorials in "),
                                    crayon_key(cli::style_hyperlink(text = "SCpubr's book",
                                                                    url = "https://enblacar.github.io/SCpubr-book/")),
                                    crayon_body(".")))

      cli::cli_text(" ")

      cli::cli_alert_success(paste0(crayon_body("If you use "),
                                    crayon_key("SCpubr"),
                                    crayon_body(" in your research, please "),
                                    crayon_key(cli::style_hyperlink(text = "cite it accordingly",
                                                                    url = "https://www.biorxiv.org/content/10.1101/2022.02.28.482303v1")),
                                    crayon_body(".")))
      cli::cli_text("  ")

      cli::cli_text(paste0(cli::style_bold(cli::col_yellow(cli::symbol$star)),
                    crayon_body(" If the package is useful to you, consider leaving a "),
                                 crayon_key("Star"),
                                 crayon_body(" in the "),
                                 crayon_key(cli::style_hyperlink(text = "GitHub repository",
                                                                 url = "https://github.com/enblacar/SCpubr")),
                                 crayon_body(".")))

      cli::cli_text(" ")

      cli::cli_text(paste0(cli::style_bold(cli::col_blue(cli::symbol$arrow_up)),
                           crayon_body(" Keep track of the package "),
                                 crayon_key("updates"),
                                 crayon_body(" on Twitter ("),
                                 crayon_key("@Enblacar"),
                                 crayon_body(") or in the "),
                                 crayon_key(cli::style_hyperlink(text = "Official NEWS website",
                                                                 url = "https://github.com/enblacar/SCpubr/blob/main/NEWS.md")),
                                 crayon_body(".")))

      cli::cli_text(" ")

      cli::cli_text(paste0(cli::style_bold(cli::col_red(cli::symbol$heart)), " ", cli::style_bold("Happy plotting!")))

      cli::cli_text(" ")

      cli::cli_h2("Checking package updates")
      cran_version <- utils::available.packages(repos = "http://cran.us.r-project.org")["SCpubr", "Version"]
      system_version <- as.character(utils::packageVersion("SCpubr"))

      l1 <- nchar(system_version)
      l2 <- nchar(cran_version)
      max_length <- max(c(l1, l2))

      message(paste0(cli::style_bold(cli::col_magenta("System: ")), cli::ansi_align(cli::col_blue(system_version), max_length, align = "right")))
      message(paste0(cli::style_bold(cli::col_magenta("CRAN:   ")), cli::ansi_align(cli::col_blue(cran_version), max_length, align = "right")))

      cli::cli_text((" "))

      if (system_version < cran_version){
         cli::cli_alert_warning(paste0(crayon_body("There is a "),
                                       crayon_key("new version"),
                                       crayon_body(" available on"),
                                       crayon_key("CRAN"),
                                       crayon_body("!")))
      } else {
        cli::cli_alert_success(paste0(crayon_body("Installation is "),
                                      crayon_key("up to date"),
                                      crayon_body(" with latest "),
                                      crayon_key("CRAN"),
                                      crayon_body(" version!")))
      }
      
      cli::cli_text(" ")

      packages <- sort(unique(unlist(SCpubr::state_dependencies(return_dependencies = TRUE))))
      functions <- sort(unique(names(SCpubr::state_dependencies(return_dependencies = TRUE))))
      functions <- sapply(functions, SCpubr:::check_suggests, passive = TRUE)
      functions <- functions[names(functions) != "Essentials"]

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
      cli::cli_text(" ")

      cli::cli_alert_danger(paste0(crayon_body("To suppress this startup message, use: "), cli::style_italic(cli::col_blue('suppressMessages(library("SCpubr"))'))))

      cli::cli_text("  ")
      message(cli::rule(col = "cyan"))
}