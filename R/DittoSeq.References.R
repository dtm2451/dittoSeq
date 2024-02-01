#' dittoSeq
#'
#' @name dittoSeq
#' @author Daniel Bunis
#' @description This package was built to make the visualization of single-cell and bulk RNA-sequencing data pipeline-agnostic and accessible for both experienced and novice coders, and for color vision impaired individuals.
#' @details Includes many plotting functions (\code{\link{dittoPlot}}, \code{\link{dittoDimPlot}}, \code{\link{dittoBarPlot}}, \code{\link{dittoHeatmap}}, ...),
#' helper funtions (\code{\link{meta}}, \code{\link{gene}}, \code{\link{isMeta}}, \code{\link{getMetas}}, ...),
#' and color adjustment functions (\code{\link{Simulate}}, \code{\link{Darken}}, \code{\link{Lighten}}),
#' to aid in making sense of RNA sequencing data.
#' All included plotting functions produce a ggplot object (or \code{\link[pheatmap]{pheatmap}} / \code{\link[ComplexHeatmap]{Heatmap}} for dittoHeatmap) by default and can spit out a full plot with just a few arguments.
#' Many additional arguments are available for customization to generate complex, publication-ready figures.
#'
#' Default \code{\link{dittoColors}} are color blindness friendly and adapted from \href{https://www.nature.com/articles/nmeth.1618}{Wong B, "Points of view: Color blindness." Nature Methods, 2011.}
#' 
#' To report bugs, suggest new features, or ask for help, the best method is to create an issue on the github, \href{https://github.com/dtm2451/dittoSeq}{here}, or the bioconductor support site (be sure to tag 'dittoSeq' so that I get a notification!), \href{https://support.bioconductor.org/}{here}
#' @importFrom utils modifyList
"_PACKAGE"

#' demuxlet.example
#'
#' @name demuxlet.example
#' @author Daniel Bunis
#' @description A dataframe containing mock demuxlet information for the 80-cell Seurat::pbmc_small dataset
#' @return A dataframe
#' @details This data was created based on the structure of real demuxlet.best output files.
#' Barcodes from Seurat's pbmc_small example data were used as the BARCODES column.
#' Cells were then assigned randomly as either SNG (singlets), DBL (doublets), or AMB (ambiguous).
#' Cells were then randomly assign to sample1-10 (or multiple samples for doublets), and this information was combined using the \code{paste} function into the typical structure of a demuxlet CALL column.
#' Random sampling of remaining data from a separate, actual, demuxlet daatset was used for remaining columns.
#' @note This is a slightly simplified example. Real demuxlet.best data has additional columns.
"demuxlet.example"
