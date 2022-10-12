% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_LigandReceptorPlot.R
\name{do_LigandReceptorPlot}
\alias{do_LigandReceptorPlot}
\title{Visualize Ligand-Receptor analysis output.}
\usage{
do_LigandReceptorPlot(
  liana_output,
  split.by = NULL,
  keep_source = NULL,
  keep_target = NULL,
  top_interactions = 25,
  dot_border = TRUE,
  border.color = "black",
  x_labels_angle = 45,
  rotate_strip_text = FALSE,
  legend.position = "bottom",
  legend.type = "colorbar",
  legend.length = 20,
  legend.width = 1,
  legend.framecolor = "grey50",
  legend.tickcolor = "white",
  legend.framewidth = 1.5,
  legend.tickwidth = 1.5,
  viridis_color_map = "G",
  viridis_direction = 1,
  font.size = 14,
  dot.size = 1,
  font.type = "sans",
  flip = FALSE,
  plot.grid = TRUE,
  grid.color = "grey90",
  grid.type = "dotted",
  compute_ChordDiagrams = FALSE,
  add_missing_LR_combinations = TRUE
)
}
\arguments{
\item{liana_output}{\strong{\code{\link[tibble]{tibble}}} | Object resulting from running \link[liana]{liana_wrap}.}

\item{split.by}{\strong{\code{\link[base]{character}}} | Whether to further facet the plot on the y axis by common ligand.complex or receptor.complex. Values to provide: NULL, ligand.complex, receptor.complex.}

\item{keep_source, keep_target}{\strong{\code{\link[base]{character}}} | Identities to keep for the source/target of the interactions. NULL otherwise.}

\item{top_interactions}{\strong{\code{\link[base]{numeric}}} | Number of unique interactions to retrieve ordered by magnitude and specificity. It does not necessarily mean that the output will contain as many, but rather an approximate value.}

\item{dot_border}{\strong{\code{\link[base]{logical}}} | Whether to draw a black border in the dots.}

\item{border.color}{\strong{\code{\link[base]{character}}} | Color to use for the border of the cells.}

\item{x_labels_angle}{\strong{\code{\link[base]{numeric}}} | One of 0 (horizontal), 45 (diagonal), 90 (vertical). Adjusts to 0 if flip = FALSE and 45 if flip = TRUE.}

\item{rotate_strip_text}{\strong{\code{\link[base]{logical}}} | Whether the text in the strips should be flipped 90 degrees.}

\item{legend.position}{\strong{\code{\link[base]{character}}} | Position of the legend in the plot. One of:
\itemize{
\item \emph{\code{top}}: Top of the figure.
\item \emph{\code{bottom}}: Bottom of the figure.
\item \emph{\code{left}}: Left of the figure.
\item \emph{\code{right}}: Right of the figure.
\item \emph{\code{none}}: No legend is displayed.
}}

\item{legend.type}{\strong{\code{\link[base]{character}}} | Type of legend to display. One of:
\itemize{
\item \emph{\code{normal}}: Default legend displayed by \pkg{ggplot2}.
\item \emph{\code{colorbar}}: Redefined colorbar legend, using \link[ggplot2]{guide_colorbar}.
\item \emph{\code{colorsteps}}: Redefined legend with colors going by range, in steps, using \link[ggplot2]{guide_colorsteps}.
}}

\item{legend.length, legend.width}{\strong{\code{\link[base]{numeric}}} | Length and width of the legend. Will adjust automatically depending on legend side.}

\item{legend.framecolor}{\strong{\code{\link[base]{character}}} | Color of the lines of the box in the legend.}

\item{legend.tickcolor}{\strong{\code{\link[base]{character}}} | Color of the ticks of the box in the legend.}

\item{legend.framewidth, legend.tickwidth}{\strong{\code{\link[base]{numeric}}} | Width of the lines of the box in the legend.}

\item{viridis_color_map}{\strong{\code{\link[base]{character}}} | A capital letter from A to H or the scale name as in \link[viridis]{scale_fill_viridis}.}

\item{viridis_direction}{\strong{\code{\link[base]{numeric}}} | Either 1 or -1. Controls how the gradient of viridis scale is formed.}

\item{font.size}{\strong{\code{\link[base]{numeric}}} | Overall font size of the plot. All plot elements will have a size relationship with this font size.}

\item{dot.size}{\strong{\code{\link[base]{numeric}}} | Size aesthetic for the dots.}

\item{font.type}{\strong{\code{\link[base]{character}}} | Base font family for the plot. One of:
\itemize{
\item \emph{\code{mono}}: Mono spaced font.
\item \emph{\code{serif}}: Serif font family.
\item \emph{\code{sans}}: Default font family.
}}

\item{flip}{\strong{\code{\link[base]{logical}}} | Whether to invert the axis using \link[ggplot2]{coord_flip}.}

\item{plot.grid}{\strong{\code{\link[base]{logical}}} | Whether to plot grid lines.}

\item{grid.color}{\strong{\code{\link[base]{character}}} | Color of the grid in the panels.}

\item{grid.type}{\strong{\code{\link[base]{character}}} | One of the possible linetype options:
\itemize{
\item \emph{\code{blank}}.
\item \emph{\code{solid}}.
\item \emph{\code{dashed}}.
\item \emph{\code{dotted}}.
\item \emph{\code{dotdash}}.
\item \emph{\code{longdash}}.
\item \emph{\code{twodash}}.
}}

\item{compute_ChordDiagrams}{\strong{\code{\link[base]{logical}}} | Whether to also compute Chord Diagrams for both the number of interactions between source and target but also between ligand.complex and receptor.complex.}

\item{add_missing_LR_combinations}{\strong{\code{\link[base]{logical}}} | Depending on the value provided to \strong{\code{top_interactions}}, there might be some source-target combinations missing. If set to TRUE, those combinations will be brought back to the plot as NA values.}
}
\value{
A ggplot2 plot with the results of the Ligand-Receptor analysis.
}
\description{
This function makes use of \href{https://github.com/saezlab/liana}{liana} package to run Ligand-Receptor analysis. Takes the output of liana and generates a dot-plot visualization according to the user's specifications.
}
\examples{
\dontrun{
  # Ligand Receptor analysis plot.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output)
  p

  # Ligand Receptor analysis plot with extra interactions.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     top_interactions = 50)
  p

  # Add missing LR combinations.
  p1 <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                      add_missing_LR_combinations = TRUE)

  # Restrict the results to only the LR combinations that have a value.
  p2 <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                      add_missing_LR_combinations = FALSE)
  p <- p1 /p2
  p

  # Ligand Receptor analysis plot with increased dot size.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     dot.size = 2)
  p

  # Ligand Receptor analysis plot with grid.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     grid.color = "steelblue")
  p

  # Ligand Receptor analysis plot with inverted axes.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     flip = TRUE,
                                     x_labels_angle = 90,
                                     rotate_strip_text = TRUE)
  p

  # Ligand Receptor analysis plot grouped by ligand.complex.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     plot.grid = TRUE,
                                     grid.color = "grey90",
                                     grid.type = "dotted",
                                     x_labels_angle = 90,
                                     rotate_strip_text = TRUE,
                                     split.by = "ligand.complex")
  p

  # Ligand Receptor analysis plot grouped by ligand.complex.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     split.by = "receptor.complex")
  p

  # Ligand Receptor analysis plot with inverted axes with filtered source and targets.
  p <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                     flip = F,
                                     keep_source = c("0", "3", "9"),
                                     keep_target = c("3", "4", "9"))
  p

  # Chord diagram of the total significant interactions from each cluster and all the rest.
  out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 25,
                                       compute_ChordDiagrams = TRUE)
  out$chord_total_interactions

  # Chord diagram of the total ligand-receptor interactions from each cluster and all the rest.
  out$chord_ligand_receptor

  # Chord diagram of the total significant interactions from each cluster and all the rest.
  out <- SCpubr::do_LigandReceptorPlot(liana_output = liana_output,
                                       top_interactions = 25,
                                       keep_source = c("0", "7"),
                                       compute_ChordDiagrams = TRUE)
  out$chord_total_interactions

  # Chord diagram of the total ligand-receptor interactions from each cluster and all the rest.
  out$chord_ligand_receptor
}
}