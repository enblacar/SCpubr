de_genes <- readRDS(system.file("extdata/de_genes_example.rds", package = "SCpubr"))
if (requireNamespace("Seurat", quietly = TRUE)) {
 suppressMessages(library("Seurat"))
}

if (requireNamespace("magrittr", quietly = TRUE)) {
 suppressMessages(library("magrittr"))
}

if (requireNamespace("dplyr", quietly = TRUE)) {
 suppressMessages(library("dplyr"))
  de_genes_scaled <- dplyr::rename(.data = de_genes,
                                   "avg_diff" = "avg_log2FC")
}

sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
metacell_mapping <- readRDS(system.file("extdata/metacell_mapping_example.rds", package = "SCpubr"))
infercnv_object <- readRDS(system.file("extdata/infercnv_object_example.rds", package = "SCpubr"))
infercnv_object_metacells <- readRDS(system.file("extdata/infercnv_object_metacells_example.rds", package = "SCpubr"))
human_chr_locations <- SCpubr::human_chr_locations
progeny_activities <- readRDS(system.file("extdata/progeny_activities_example.rds", package = "SCpubr"))
dorothea_activities <- readRDS(system.file("extdata/dorothea_activities_example.rds", package = "SCpubr"))
enriched_terms <- readRDS(system.file("extdata/enriched_terms_example.rds", package = "SCpubr"))


# Get packages.
dependencies <- SCpubr::state_dependencies(return_dependencies = TRUE)

dependencies[["utils"]] <- c("Seurat",
                             "rlang",
                             "dplyr",
                             "magrittr",
                             "dplyr",
                             "tidyr",
                             "tibble",
                             "stringr",
                             "plyr",
                             "grDevices",
                             "stats",
                             "grid",
                             "assertthat",
                             "ComplexHeatmap")

# Check them.
dep_check <- list()
for (func in names(dependencies)){
  packages <- dependencies[[func]]
  value = FALSE
  for (pkg in packages){
    if (!requireNamespace(pkg, quietly = TRUE)) {
        value <- TRUE
    }
  }
  dep_check[[func]] <- value
}


# Remove this for publication in CRAN.
if (isFALSE(dep_check[["do_LigandReceptorPlot"]])){
  liana_output <- readRDS(system.file("extdata/liana_output_example.rds", package = "SCpubr"))
}

if (isFALSE(dep_check[["do_DimPlot"]]) &
    isFALSE(dep_check[["do_CorrelationPlot"]]) &
    isFALSE(dep_check[["do_ChordDiagramPlot"]]) &
    isTRUE(requireNamespace(pkg, quietly = TRUE)) &
    isFALSE(dep_check[["save_Plot"]])){
  p <- SCpubr::do_DimPlot(sample)
  p.heatmap <- SCpubr::do_CorrelationPlot(sample)
  data <- p.heatmap@ht_list$`Pearson coef.`@matrix
  p.pheatmap <- pheatmap::pheatmap(data)
  p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
  figure_path <- getwd()
}


#monocle_sample <- sample
#monocle_cds <- test.data$monocle_cds
