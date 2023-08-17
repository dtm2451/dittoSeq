#' Returns the names of all dimensionality reduction slots of a target object.
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @return A string vector of the names of all dimensionality reduction slots of the \code{object}.
#' These represent the options for the \code{reduction.use} input of \code{\link{dittoDimPlot}}.
#'
#' @examples
#'
#' example("addDimReduction", echo = FALSE)
#'
#' # To see all metadata slots of an object
#' getReductions(myRNA)
#'
#' @author Daniel Bunis
#' @importFrom SingleCellExperiment reducedDimNames
#' @export

getReductions <- function(object){

    reds <- NULL
    if (is(object,"SingleCellExperiment")) {
        reds <- SingleCellExperiment::reducedDimNames(object)
    }
    if (is(object,"Seurat")) {
        .error_if_no_Seurat()
        reds <- Seurat::Reductions(object)
    }
    if (is(object,"seurat")) {
        reds <- names(object@dr)
    }
    
    # Standardize non-existent reductions output
    if (identical(reds, NA) || length(reds)==0) {
        reds <- NULL
    }
    
    reds
}
