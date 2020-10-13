#' Retrieve whether a SingleCellExperiment object would be treated as bulk versus single-cell by dittoSeq
#' @param object A target SingleCellExperiment object
#'
#' Alternatively, anything else, but then the result will always be \code{FALSE}
#' @return Logical: whether the provided object would be treated as bulk data by dittoSeq
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' isBulk(myRNA)
#' 
#' scRNA <- setBulk(myRNA, FALSE)
#' isBulk(scRNA)
#' @importFrom SingleCellExperiment "int_metadata<-" int_metadata
#' @export
isBulk <- function(object) {
    OUT <- FALSE
    if (is(object,"SingleCellExperiment")) {
        if (!is.null(SingleCellExperiment::int_metadata(object)$bulk)) {
            if (SingleCellExperiment::int_metadata(object)$bulk) {
                OUT <- TRUE
            }
        }
    } else if (is(object,"SummarizedExperiment")) {
        OUT <- TRUE
    }
    OUT
}

#' Set whether a SingleCellExperiment object should be treated as bulk versus single-cell by dittoSeq
#' @rdname setBulk
#' @param object A target SingleCellExperiment object
#' @param set Logical, whether the object should be considered as bulk (TRUE) or not (FALSE)
#' @return A \code{\linkS4class{SingleCellExperiment}} object with "bulk" internal metadata set to \code{set}
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' isBulk(myRNA)
#'
#' scRNA <- setBulk(myRNA, FALSE)
#' isBulk(scRNA)
#'
#' # Now, if we make a heatmap with this data, we will see that single-cell
#' # defaults (ordering by the first 'annot.by' & cell names not shown) are used.
#' dittoHeatmap(scRNA, getGenes(scRNA)[1:30],
#'     annot.by = c("clustering", "groups"),
#'     main = "isBulk(object) == FALSE")
#'
#' @importFrom SingleCellExperiment "int_metadata<-" int_metadata
#' @export
setGeneric("setBulk", function(object, set = TRUE) standardGeneric("setBulk"))

#' @rdname setBulk
setMethod("setBulk", "SingleCellExperiment", function(
    object, set = TRUE) {

    int_meta <- SingleCellExperiment::int_metadata(object)
    int_meta$bulk <- set
    SingleCellExperiment::int_metadata(object) <- int_meta

    object
})
