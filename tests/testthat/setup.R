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

# Remove this for publication in CRAN.
liana_output <- test.data$liana_output
p <- SCpubr::do_DimPlot(sample)
p.heatmap <- SCpubr::do_CorrelationPlot(sample)
data <- p.heatmap@ht_list$`Pearson coef.`@matrix
p.pheatmap <- pheatmap::pheatmap(data)
p.chord <- SCpubr::do_ChordDiagramPlot(sample = sample, from = "seurat_clusters", to = "orig.ident")
figure_path <- getwd()
monocle_sample <- sample
monocle_cds <- test.data$monocle_cds
