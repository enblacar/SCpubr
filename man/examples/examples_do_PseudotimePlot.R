\donttest{
  # Genes to use.
  pseudotime_genes <- c("CD14", "LYN")

  # Define your sample.
  # sample <- your_seurat_object
  # Transform into CDS.
  # cds <- SeuratWrappers::as.cell_data_set(sample)

  # Compute monocle clusters and partitions.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = TRUE,
                                   pseudotime_genes = pseudotime_genes)

  # Compute monocle clusters and keep a single partition.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = TRUE,
                                   pseudotime_genes = pseudotime_genes)

  # Compute monocle partitions but keep original identities as clusters.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes)

  # Keep original identities as clusters and a single partition.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes)

  # Set a metadata varible as clusters and a single partition.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = FALSE,
                                   group.by = "orig.ident")


  # Compute trajectory graph.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes)
  # Retrieve trajectory groups.
  p1 <- out$trajectory_groups
  # Retrieve trajectory partitions.
  p2 <- out$trajectory_partitions

  p <- p1 | p2
  p

  # Change trajectory graph width.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 1)
  p1 <- out$trajectory_partitions

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 2)
  p2 <- out$trajectory_partitions

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 3)
  p3 <- out$trajectory_partitions


  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 3,
                                   trajectory_graph_color = "white")
  p4 <- out$trajectory_partitions

  p <- (p1 | p2) / (p3 | p4)
  p

  # Add nodes, branches and leaves to the trajectory graph.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 2,
                                   label_roots = TRUE)
  p1 <- out$trajectory_partitions

  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 2,
                                   label_roots = TRUE,
                                   label_branches = TRUE)
  p2 <- out$trajectory_partitions


  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   trajectory_graph_segment_size = 2,
                                   label_roots = TRUE,
                                   label_branches = TRUE,
                                   label_leaves = TRUE)
  p3 <- out$trajectory_partitions

  p <- p1 | p2 | p3
  p

  # Plot pseudotime with monocle partitions using highest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   is_max_score_the_start = TRUE)
  p1 <- out$pseudotime

  # Plot pseudotime with monocle partitions using lowest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   is_max_score_the_start = FALSE)
  p2 <- out$pseudotime


  # Plot pseudotime with monocle partitions using highest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   is_max_score_the_start = TRUE)
  p3 <- out$pseudotime

  # Plot pseudotime with monocle partitions using lowest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = FALSE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   is_max_score_the_start = FALSE)
  p4 <- out$pseudotime

  p <- (p1 | p2) / (p3 | p4)
  p

  # Plot pseudotime with monocle partitions using highest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   symmetrical_scale = FALSE,
                                   is_max_score_the_start = TRUE)
  p1 <- out$pseudotime
  p3 <- out$enrichment
  # Plot pseudotime with monocle partitions using lowest score as root.
  out <- SCpubr::do_PseudotimePlot(sample = sample,
                                   cds = cds,
                                   compute_monocle_partitions = TRUE,
                                   compute_monocle_clusters = FALSE,
                                   pseudotime_genes = pseudotime_genes,
                                   symmetrical_scale = TRUE,
                                   is_max_score_the_start = FALSE)
  p2 <- out$pseudotime
  p4 <- out$enrichment

  p <- (p1 | p2) / (p3 | p4)
  p
}
