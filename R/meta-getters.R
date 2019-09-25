#### is.meta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot in a target object.
#'
#' @param test String or vector of strings, the "potential.metadata.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a meta.data slot within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed metadata slots if \code{return.values = TRUE}.
#' @details
#' For Seurat objects, also returns TRUE for the input \code{"ident"} because, for all dittoSeq visualiztions, \code{"ident"} will retrieve a Seurat objects' clustering slot.
#'
#' @seealso
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all metadata slots of an object
#' get.metas(pbmc)
#'
#' # To check if something is a metadata slot
#' is.meta("age", object = pbmc) # False
#' is.meta("nCount_RNA", object = pbmc) # True
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' is.meta("groups")
#'
#' # works for multiple test metas
#' is.meta(c("age","nCount_RNA"))
#'
#' # and with return.values = TRUE, returns all elements that are indeed slots
#' is.meta(c("age","nCount_RNA", "RNA_snn_res.0.8"),
#'     return.values = TRUE)
#'
#' @export
#' @import ggplot2
#' @importFrom utils packageVersion

is.meta <- function(test, object=DEFAULT, return.values=FALSE){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    if (return.values) {
        return(test[is.meta(test, object, return.values=FALSE)])
    } else {
        metas <- get.metas(object)
        if (grepl("Seurat", .class_of(object))) {
            metas <- c(metas,"ident")
        }
        return(test %in% metas)
    }
}

#### get.metas: prints the names of all the metadata lists for the object ####
#' Returns the names of all meta.data slots of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @param names.only Logical, \code{TRUE} by default, which sets whether just the names should be output versus the entire metadata dataframe.
#' @return A string vector of the names of all metadata slots of the \code{object}, or alternatively the entire dataframe of metadatas if \code{names.only} is set to \code{FALSE}
#' @seealso
#' \code{\link{is.meta}} for checking if certain metadata slots exist in an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all metadata slots of an object
#' get.metas(pbmc)
#'
#' # To retrieve the entire metadata matrix
#' get.metas(pbmc, names.only = FALSE)
#'
#' @export

get.metas <- function(object=DEFAULT, names.only = TRUE){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    metadata <-
        if (.class_of(object)=="SingleCellExperiment") {
            SummarizedExperiment::colData(object)
        } else {
            object@meta.data
        }
    if (names.only) {
        return(names(metadata))
    } else {
        return(metadata)
    }
}

#### meta: for extracting the values of a particular metadata for all cells/samples ####
#' Returns the values of a meta.data  for all cells/samples
#'
#' @param meta String, the name of the "metadata" slot to grab. OR "ident" to retireve the clustering of a Seurat \code{object}.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @return Returns the values of a metadata slot, or the clustering slot if \code{meta = "ident"} and the \code{object} is a Seurat.
#' @seealso
#' \code{\link{meta.levels}} for returning just the unique discrete identities that exist within a metadata slot
#'
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{is.meta}} for testing whether something is the name of a metadata slot
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#' meta("RNA_snn_res.1", object = "pbmc")
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta("RNA_snn_res.1")
#'
#' @export

meta <- function(meta, object=DEFAULT){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }

    if( meta=="ident" & grepl("Seurat", .class_of(object))) {
    # Retrieve clustering from Seurats
        if (grepl("v3",.class_of(object))) {
            return(as.character(Seurat::Idents(object)))
        } else {
            return(as.character(
                eval(expr = parse(text = paste0(
                    "object@ident")))))
        }
    }
    get.metas(object, names.only = FALSE)[,meta]
}

#### meta.levels: for obtaining the different classifications of a meta.data
#' Gives the distinct values of a meta.data slot (or ident)
#'
#' @param meta quoted "meta.data.slot" name = REQUIRED. the meta.data slot whose potential values should be retrieved.
#' @param object the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @return Returns the distinct values of a metadata slot given to all cells/samples or for a subset of cells/samples.
#' (Alternatively, returns the distinct values of clustering if \code{meta = "ident"} and the object is a \code{Seurat} object).
#' @seealso
#' \code{\link{meta}} for returning an entire metadata slots of an \code{object}, not just the potential levels
#'
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{is.meta}} for testing whether something is the name of a metadata slot
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' meta.levels("RNA_snn_res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta.levels("RNA_snn_res.1")
#' @export

meta.levels <- function(meta, object = DEFAULT, cells.use = NULL){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    meta.values <- as.character(meta(meta, object))
    if (!is.null(cells.use)) {
        all.cells <- .all_cells(object)
        cells.use <- .which_cells(cells.use, object)
        meta.values <- meta.values[all.cells %in% cells.use]
    }
    levels(as.factor(meta.values))
}
