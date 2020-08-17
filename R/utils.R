#' @importFrom SingleCellExperiment SingleCellExperiment
.make_sce_if_raw <- function(object, metadata = NULL, reduction.matrix = NULL) {
    
    # Return if object is valid already.
    if (is(object, "SingleCellExperiment") ||
        is(object, "Seurat") ||
        is(object, "seurat")) {
        
        return(object)
    }
    
    meta_given <- !is.null(metadata)
    reduction_given <- !is.null(reduction.matrix)
    
    ### Validate colnames
    # Use cell names of object if there & unique, else from rownames of metadata.
    if (is.null(colnames) || ncol(object) != length(unique(colnames(object)))) {
        part1 <- "colnames('object') are NULL or invalid,"
        if (meta_given) {
            message(part1, " using rownames('metadata').")
            colnames(object) <- rownames(metadata)
        } else {
            message(part1, " using col1, col2, col3, etc.")
            colnames(object) <- paste0("col", seq_len(ncol(object)))
        }
    } else {
        if (meta_given) {
            rownames(metadata) <- colnames(object)
        }
    }
    
    ### Make SCE
    args <- list(assays = list(counts = object))
    if (meta_given) {
        args$colData <- metadata
    }
    if (reduction_given) {
        args$reducedDims <- list(Dim = reduction.matrix)
    }
    
    do.call(SingleCellExperiment::SingleCellExperiment, args)
}

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
