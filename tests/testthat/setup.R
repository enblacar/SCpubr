suppressMessages(library(Seurat))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
test.data <- SCpubr:::test.data
sample <- test.data$sample
metacell_mapping <- test.data$metacell_mapping
infercnv_object <- test.data$infercnv_object
infercnv_object_metacells <- test.data$infercnv_object_metacells
human_chr_locations <- test.data$human_chr_locations
de_genes <- test.data$de_genes
de_genes_scaled <- de_genes %>% dplyr::rename("avg_diff" = "avg_log2FC")
progeny_activities <- test.data$progeny_activities
dorothea_activities <- test.data$dorothea_activities
enriched_terms <- test.data$enriched_terms
