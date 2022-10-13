suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))

sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))
metacell_mapping <- readRDS(system.file("extdata/metacell_mapping_example.rds", package = "SCpubr"))
infercnv_object <- readRDS(system.file("extdata/infercnv_object_example.rds", package = "SCpubr"))
infercnv_object_metacells <- readRDS(system.file("extdata/infercnv_object_metacells_example.rds", package = "SCpubr"))
human_chr_locations <- SCpubr::human_chr_locations
de_genes <- readRDS(system.file("extdata/de_genes_example.rds", package = "SCpubr"))
de_genes_scaled <- de_genes %>% dplyr::rename("avg_diff" = "avg_log2FC")
progeny_activities <- readRDS(system.file("extdata/progeny_activities_example.rds", package = "SCpubr"))
dorothea_activities <- readRDS(system.file("extdata/dorothea_activities_example.rds", package = "SCpubr"))
enriched_terms <- readRDS(system.file("extdata/enriched_terms_example.rds", package = "SCpubr"))




# Remove this for publication in CRAN.
liana_output <- readRDS(system.file("extdata/liana_output_example.rds", package = "SCpubr"))
p <- SCpubr::do_DimPlot(sample)
p.heatmap <- SCpubr::do_CorrelationPlot(sample)
data <- p.heatmap@ht_list$`Pearson coef.`@matrix
p.pheatmap <- pheatmap::pheatmap(data)
p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
figure_path <- getwd()
# monocle_sample <- sample
# monocle_cds <- test.data$monocle_cds
