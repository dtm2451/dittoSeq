#' @importFrom SingleCellExperiment reducedDim int_metadata
#' @importFrom SummarizedExperiment assay assays
#' @importFrom cowplot ggdraw get_legend

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
