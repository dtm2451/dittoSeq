#' DittoSeq
#'
#' @docType package
#' @name DittoSeq
#' @description This package was built to make the analysis and visualization of single-cell and bulk RNA-sequencing data accessible for both experience and novice coders, and for colorblind individuals.
#' @details Includes many plotting functions (DBPlot, DBDimPlot, DBBarPlot, ...), color adjustment functions (Simulate, Darken, Lighten), and helper funtions (meta, gene, is.meta, get.metas,...) to aid in making sense of single cell or bulk RNA sequencing data.  All included plotting functions produce a ggplot and can spit out full plot with just a few arguments.  Many arguments are available for customization to generate complex publication-ready figures.  For more information, to give feedback, or to suggest new features, see the github, [here](https://github.com/dtm2451/DittoSeq). Default color panel, stored in MYcolors, is colorblind friendly [Wong B, "Points of view: Color blindness." Nature Methods, 2011.](https://www.nature.com/articles/nmeth.1618).
NULL

#' Protanope and Deuteranope friendly colors adapted from Wong, B. Nature Methods, 2011.
#' @author Daniel Bunis
#' @description A Character vector of 24 colors in hexadecimal form. Default color.panel for DittoSeq Visualizations
#'
#' 1-7 = Suggested color panel from Wong, B. Nature Methods, 2011;
#' 8 = gray40;
#' 9-16 = 25% darker versions of colors 1-8;
#' 17-24 = 25% lighter versions of colors 1-8
#'
"MYcolors"


