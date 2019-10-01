#' dittoSeq
#'
#' @docType package
#' @name dittoSeq
#' @description This package was built to make the analysis and visualization of single-cell and bulk RNA-sequencing data accessible for both experience and novice coders, and for colorblind individuals.
#' @details Includes many plotting functions (\code{\link{dittoPlot}}, \code{\link{dittoDimPlot}}, \code{\link{dittoBarPlot}}, \code{\link{dittoHeatmap}}, ...),
#' color adjustment functions (\code{\link{Simulate}}, \code{\link{Darken}}, \code{\link{Lighten}}),
#' and helper funtions (\code{\link{meta}}, \code{\link{gene}}, \code{\link{is.meta}}, \code{\link{get.metas}}, ...) to aid in making sense of single cell or bulk RNA sequencing data.
#' All included plotting functions produce a ggplot (or pheatmap for dittoHeatmap) and can spit out full plot with just a few arguments.
#' Many additional arguments are available for customization to generate complex publication-ready figures.
#'
#' Default color panel is colorblind friendly [Wong B, "Points of view: Color blindness." Nature Methods, 2011.](https://www.nature.com/articles/nmeth.1618).
#'
#' For more information, to give feedback, or to suggest new features, see the github, [here](https://github.com/dtm2451/DittoSeq).
NULL
