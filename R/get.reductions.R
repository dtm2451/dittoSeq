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

    if (is(object,"SingleCellExperiment")) {
        return(SingleCellExperiment::reducedDimNames(object))
    }
    if (is(object,"Seurat")) {
        .error_if_no_Seurat()
        return(Seurat::Reductions(object))
    }
    if (is(object,"seurat")) {
        return(names(object@dr))
    }
    
    NULL
}
