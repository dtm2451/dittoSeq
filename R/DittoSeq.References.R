#' dittoSeq
#'
#' @docType package
#' @name dittoSeq
#' @author Daniel Bunis
#' @description This package was built to make the analysis and visualization of single-cell and bulk RNA-sequencing data accessible for both experience and novice coders, and for colorblind individuals.
#' @details Includes many plotting functions (\code{\link{dittoPlot}}, \code{\link{dittoDimPlot}}, \code{\link{dittoBarPlot}}, \code{\link{dittoHeatmap}}, ...),
#' color adjustment functions (\code{\link{Simulate}}, \code{\link{Darken}}, \code{\link{Lighten}}),
#' and helper funtions (\code{\link{meta}}, \code{\link{gene}}, \code{\link{isMeta}}, \code{\link{getMetas}}, ...) to aid in making sense of single cell or bulk RNA sequencing data.
#' All included plotting functions produce a ggplot (or plotly, or pheatmap for dittoHeatmap) and can spit out full plot with just a few arguments.
#' Many additional arguments are available for customization to generate complex publication-ready figures.
#'
#' Default color panel is colorblind friendly [Wong B, "Points of view: Color blindness." Nature Methods, 2011.](https://www.nature.com/articles/nmeth.1618).
#'
#' For more information, to give feedback, or to suggest new features, see the github, [here](https://github.com/dtm2451/DittoSeq).
NULL

#' demuxlet.example
#'
#' @name demuxlet.example
#' @author Daniel Bunis
#' @description A dataframe containing mock demuxlet information for the 80-cell Seurat::pbmc_small dataset
#' @return A dataframe
#' @details This data was created based on the structure of real demuxlet.best output files.
#' Barcodes from the \code{\link[Seurat]{pbmc_small}} dataset were used as the BARCODES column.
#' Cells were then assigned randomly as either SNG (singlets), DBL (doublets), or AMB (ambiguous).
#' Cells were then randomly assign to sample1-10 (or multiple samples for doublets), and this information was combined using the \code{paste} function into the typical structure of a demuxlet CALL column.
#' Random sampling of remaining data from a separate, actual, demuxlet daatset was used for remaining columns.
#' @note This is a slightly simplified example. Real demuxlet.best data has additional columns.
"demuxlet.example"
