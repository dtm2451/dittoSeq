.msg_if <- function(verbose, ...){
    if (verbose) {
        message(...)
    }
}

.error_if_no_Seurat <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat installation required for working with Seurat objects")
    }
}

.error_if_no_plotly <- function() {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("plotly installation required for using hover")
    }
}

.error_if_no_ggplot.multistats <- function() {
    if (!requireNamespace("ggplot.multistats", quietly = TRUE)) {
        stop("ggplot.multistats installation required for supplying 'color.var' to dittoHex plotters.")
    }
}

.error_if_no_ggrastr <- function() {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
        stop("ggrastr installation required for using rasterization with dittoScatterPlot plotters.")
    }
}

.error_if_no_complexHm <- function() {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("ComplexHeatmap installation required for using `complex` in dittoHeatmap.",
             "\nInstall with BiocManager::install('ComplexHeatmap').")
    }
}
