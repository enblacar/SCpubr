% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_CellularStatesPlot.R
\name{do_CellularStatesPlot}
\alias{do_CellularStatesPlot}
\title{Cellular States plot.}
\usage{
do_CellularStatesPlot(
  sample,
  input_gene_list,
  x1,
  y1,
  x2 = NULL,
  y2 = NULL,
  group.by = NULL,
  colors.use = NULL,
  colorblind = FALSE,
  legend.position = "bottom",
  legend.icon.size = 4,
  legend.ncol = NULL,
  legend.nrow = NULL,
  legend.byrow = FALSE,
  plot.title = NULL,
  plot.subtitle = NULL,
  plot.caption = NULL,
  font.size = 14,
  font.type = "sans",
  xlab = NULL,
  ylab = NULL,
  axis.ticks = TRUE,
  axis.text = TRUE,
  verbose = FALSE,
  enforce_symmetry = FALSE,
  plot_marginal_distributions = FALSE,
  marginal.type = "density",
  marginal.size = 5,
  marginal.group = TRUE,
  plot_cell_borders = TRUE,
  plot_enrichment_scores = FALSE,
  border.size = 2,
  border.color = "black",
  pt.size = 2,
  raster = FALSE,
  raster.dpi = 1024,
  plot_features = FALSE,
  features = NULL,
  use_viridis = TRUE,
  viridis.palette = "G",
  viridis.direction = 1,
  sequential.palette = "YlGnBu",
  sequential.direction = -1,
  nbin = 24,
  ctrl = 100,
  number.breaks = 5,
  plot.title.face = "bold",
  plot.subtitle.face = "plain",
  plot.caption.face = "italic",
  axis.title.face = "bold",
  axis.text.face = "plain",
  legend.title.face = "bold",
  legend.text.face = "plain"
)
}
\arguments{
\item{sample}{\strong{\code{\link[SeuratObject]{Seurat}}} | A Seurat object, generated by \link[Seurat]{CreateSeuratObject}.}

\item{input_gene_list}{\strong{\code{\link[SCpubr]{named_list}}} | Named list of lists of genes to be used as input.}

\item{x1}{\strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. First feature in the X axis. Will go on the right side of the X axis if y2 is not provided and top-right quadrant if provided.}

\item{y1}{\strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. First feature on the Y axis. Will become the Y axis if y2 is not provided and bottom-right quadrant if provided.}

\item{x2}{\strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. Second feature on the X axis. Will go on the left side of the X axis if y2 is not provided and top-left quadrant if provided.}

\item{y2}{\strong{\code{\link[base]{character}}} | A name of a list from input_gene_list. Second feature on the Y axis. Will become the bottom-left quadrant if provided.}

\item{group.by}{\strong{\code{\link[base]{character}}} | Metadata variable to group the output by. Has to be a character of factor column.}

\item{colors.use}{\strong{\code{\link[SCpubr]{named_vector}}} | Named vector of valid color representations (either name of HEX codes) with as many named colors as unique values of group.by. If group.by is not provided, defaults to the unique values of \link[Seurat]{Idents}. If not provided, a color scale will be set by default.}

\item{colorblind}{\strong{\code{\link[base]{logical}}} | Whether to use colorblind-friendly colors for categorical variables. In place when \code{colors.use} is not used. Allows for a maximum of 85 different classes within a categorical variable.}

\item{legend.position}{\strong{\code{\link[base]{character}}} | Position of the legend in the plot. One of:
\itemize{
\item \emph{\code{top}}: Top of the figure.
\item \emph{\code{bottom}}: Bottom of the figure.
\item \emph{\code{left}}: Left of the figure.
\item \emph{\code{right}}: Right of the figure.
\item \emph{\code{none}}: No legend is displayed.
}}

\item{legend.icon.size}{\strong{\code{\link[base]{numeric}}} | Size of the icons in legend.}

\item{legend.ncol}{\strong{\code{\link[base]{numeric}}} | Number of columns in the legend.}

\item{legend.nrow}{\strong{\code{\link[base]{numeric}}} | Number of rows in the legend.}

\item{legend.byrow}{\strong{\code{\link[base]{logical}}} | Whether the legend is filled by row or not.}

\item{plot.title, plot.subtitle, plot.caption}{\strong{\code{\link[base]{character}}} | Title, subtitle or caption to use in the plot.}

\item{font.size}{\strong{\code{\link[base]{numeric}}} | Overall font size of the plot. All plot elements will have a size relationship with this font size.}

\item{font.type}{\strong{\code{\link[base]{character}}} | Base font family for the plot. One of:
\itemize{
\item \emph{\code{mono}}: Mono spaced font.
\item \emph{\code{serif}}: Serif font family.
\item \emph{\code{sans}}: Default font family.
}}

\item{xlab, ylab}{\strong{\code{\link[base]{character}}} | Titles for the X and Y axis.}

\item{axis.ticks}{\strong{\code{\link[base]{logical}}} | Whether to show axis ticks.}

\item{axis.text}{\strong{\code{\link[base]{logical}}} | Whether to show axis text.}

\item{verbose}{\strong{\code{\link[base]{logical}}} | Whether to show extra comments, warnings,etc.}

\item{enforce_symmetry}{\strong{\code{\link[base]{logical}}} | Whether to enforce the plot to follow a symmetry (3 variables, the X axis has 0 as center, 4 variables, all axis have the same range and the plot is squared).}

\item{plot_marginal_distributions}{\strong{\code{\link[base]{logical}}} |  Whether to plot marginal distributions on the figure or not.}

\item{marginal.type}{\strong{\code{\link[base]{character}}} | One of:
\itemize{
\item \emph{\code{density}}: Compute density plots on the margins.
\item \emph{\code{histogram}}: Compute histograms on the margins.
\item \emph{\code{boxplot}}: Compute boxplot on the margins.
\item \emph{\code{violin}}: Compute violin plots on the margins.
\item \emph{\code{densigram}}: Compute densigram plots on the margins.
}}

\item{marginal.size}{\strong{\code{\link[base]{numeric}}} | Size ratio between the main and marginal plots. A value of 5 means that the main plot is 5 times bigger than the marginal plots.}

\item{marginal.group}{\strong{\code{\link[base]{logical}}} | Whether to group the marginal distribution by group.by or current identities.}

\item{plot_cell_borders}{\strong{\code{\link[base]{logical}}} | Whether to plot border around cells.}

\item{plot_enrichment_scores}{\strong{\code{\link[base]{logical}}} | Whether to report enrichment scores for the input lists as plots.}

\item{border.size}{\strong{\code{\link[base]{numeric}}} | Width of the border of the cells.}

\item{border.color}{\strong{\code{\link[base]{character}}} | Color for the border of the heatmap body.}

\item{pt.size}{\strong{\code{\link[base]{numeric}}} | Size of the dots.}

\item{raster}{\strong{\code{\link[base]{logical}}} | Whether to raster the resulting plot. This is recommendable if plotting a lot of cells.}

\item{raster.dpi}{\strong{\code{\link[base]{numeric}}} | Pixel resolution for rasterized plots. Defaults to 1024. Only activates on Seurat versions higher or equal than 4.1.0.}

\item{plot_features}{\strong{\code{\link[base]{logical}}} | Whether to also report any other feature onto the primary plot.}

\item{features}{\strong{\code{\link[base]{character}}} | Additional features to plot.}

\item{use_viridis}{\strong{\code{\link[base]{logical}}} | Whether to use viridis color scales.}

\item{viridis.palette}{\strong{\code{\link[base]{character}}} | A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.}

\item{viridis.direction}{\strong{\code{\link[base]{numeric}}} | Either 1 or -1. Controls how the gradient of viridis scale is formed.}

\item{sequential.palette}{\strong{\code{\link[base]{character}}} | Type of sequential color palette to use. Out of the sequential palettes defined in \strong{\code{\link[RColorBrewer]{brewer.pal}}}.}

\item{sequential.direction}{\strong{\code{\link[base]{numeric}}} | Direction of the sequential color scale. Either 1 or -1.}

\item{nbin}{\strong{\code{\link[base]{numeric}}} | Number of bins to use in \link[Seurat]{AddModuleScore}.}

\item{ctrl}{\strong{\code{\link[base]{numeric}}} | Number of genes in the control set to use in \link[Seurat]{AddModuleScore}.}

\item{number.breaks}{\strong{\code{\link[base]{numeric}}} | Controls the number of breaks in continuous color scales of ggplot2-based plots.}

\item{plot.title.face, plot.subtitle.face, plot.caption.face, axis.title.face, axis.text.face, legend.title.face, legend.text.face}{\strong{\code{\link[base]{character}}} | Controls the style of the font for the corresponding theme element.  One of:
\itemize{
\item \emph{\code{plain}}: For normal text.
\item \emph{\code{italic}}: For text in itallic.
\item \emph{\code{bold}}: For text in bold.
\item \emph{\code{bold.italic}}: For text both in itallic and bold.
}}
}
\value{
A ggplot2 object containing a butterfly plot.
}
\description{
This plot aims to show the relationships between distinct enrichment scores. If 3 variables are provided, the relationship is between the Y axis and the dual X axis.
If 4 variables are provided, each corner of the plot represents how enriched the cells are in that given list. How to interpret this? In a 3-variable plot, the Y axis
just means one variable. The higher the cells are in the Y axis the more enriched they are in that given variable. The X axis is a dual parameter one. Cells falling
into each extreme of the axis are highly enriched for either x1 or x2, while cells falling in between are not enriched for any of the two. In a 4-variable plot, each corner
shows the enrichment for one of the 4 given features. Cells will tend to locate in either of the four corners, but there will be cases of cells locating mid-way between two
given corners (enriched in both features) or in the middle of the plot (not enriched for any).
}
\details{
This plots are based on the following publications:
\itemize{
\item Neftel, C. \emph{et al}. An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma. Cell 178, 835-849.e21 (2019). \doi{10.1016/j.cell.2019.06.024}
\item Tirosh, I., Venteicher, A., Hebert, C. \emph{et al}. Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma. Nature 539, 309–313 (2016). \doi{10.1038/nature20123}
}
}
\examples{
\donttest{
  # Check Suggests.
  value <- SCpubr:::check_suggests(function_name = "do_CellularStatesPlot", passive = TRUE)

  if (isTRUE(value)){
    # Consult the full documentation in https://enblacar.github.io/SCpubr-book/

    # Define your Seurat object.
    sample <- readRDS(system.file("extdata/seurat_dataset_example.rds", package = "SCpubr"))

    # Define some gene sets to query. It has to be a named list.
    gene_set <- list("A" = rownames(sample)[1:10],
                     "B" = rownames(sample)[11:20],
                     "C" = rownames(sample)[21:30],
                     "D" = rownames(sample)[31:40])

    # Using two variables: A scatter plot X vs Y.
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "B",
                                       nbin = 1,
                                       ctrl = 10)
    p


    # Using three variables. Figure from: https://www.nature.com/articles/nature20123.
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "B",
                                       x2 = "C",
                                       nbin = 1,
                                       ctrl = 10)
    p


    # Using four variables. Figure from: https://pubmed.ncbi.nlm.nih.gov/31327527/
    p <- SCpubr::do_CellularStatesPlot(sample = sample,
                                       input_gene_list = gene_set,
                                       x1 = "A",
                                       y1 = "C",
                                       x2 = "B",
                                       y2 = "D",
                                       nbin = 1,
                                       ctrl = 10)
    p
  } else if (base::isFALSE(value)){
    message("This function can not be used without its suggested packages.")
    message("Check out which ones are needed using `SCpubr::state_dependencies()`.")
  }
}

}
