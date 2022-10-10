suppressMessages(library(Seurat))
suppressMessages(library(monocle3))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
test.data <- SCpubr:::test.data
sample <- test.data$sample
metacell_mapping <- test.data$metacell_mapping
infercnv_object <- test.data$infercnv_object
infercnv_object_metacells <- test.data$infercnv_object_metacells
human_chr_locations <- test.data$human_chr_locations
de_genes <- test.data$de_genes
de_genes_scaled <- test.data$de_genes_scaled
liana_output <- test.data$liana_output
progeny_activities <- test.data$progeny_activities
dorothea_activities <- test.data$dorothea_activities
monocle_sample <- test.data$monocle_sample
monocle_cds <- test.data$monocle_cds
enriched_terms <- test.data$enriched_terms
