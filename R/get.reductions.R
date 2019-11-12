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

    if (.class_of(object)=="SingleCellExperiment") {
        reductions <- SingleCellExperiment::reducedDimNames(object)
    }
    if (.class_of(object)=="Seurat.v3") {
        reductions <- Seurat::Reductions(object)
    }
    if (.class_of(object)=="Seurat.v2") {
        reductions <- names(object@dr)
    }
    if (.class_of(object)=="RNAseq") {
        reductions <- names(object@reductions)
    }

    reductions
}
