suppressMessages(library(Seurat))
suppressMessages(library(infercnv))
suppressMessages(library(monocle3))
suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
test.data <- readRDS(file= system.file("test_data/test.data.rds", package= "SCpubr"))
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
p <- SCpubr::do_DimPlot(sample)
p.heatmap <- SCpubr::do_CorrelationPlot(sample)
data <- p.heatmap@ht_list$`Pearson coef.`@matrix
p.pheatmap <- pheatmap::pheatmap(data)
p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
figure_path <- getwd()
