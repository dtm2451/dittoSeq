#' Returns the names of all dimensionality reduction slots of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @return A string vector of the names of all dimensionality reduction slots of the \code{object}.
#' These represent the options for the \code{reduction.use} input of \code{\link{dittoDimPlot}}.
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all metadata slots of an object
#' getReductions(pbmc)
#'
#' @author Daniel Bunis
#' @export

getReductions <- function(object=DEFAULT){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }

    if (is(object,"SingleCellExperiment")) {
        return(SingleCellExperiment::reducedDimNames(object))
    }
    if (is(object,"Seurat")) {
        return(Seurat::Reductions(object))
    }
    if (is(object,"seurat")) {
        return(names(object@dr))
    }
}
