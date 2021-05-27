#### isMeta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot in a target object.
#'
#' @param test String or vector of strings, the "potential.metadata.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a meta.data slot within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed metadata slots if \code{return.values = TRUE}.
#' @details
#' For Seurat objects, also returns TRUE for the input \code{"ident"} because, for all dittoSeq visualiztions, \code{"ident"} will retrieve a Seurat objects' clustering slot.
#'
#' @seealso
#' \code{\link{getMetas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#'
#' # To check if something is a metadata slot
#' isMeta("timepoint", object = myRNA) # FTRUE
#' isMeta("nCount_RNA", object = myRNA) # FALSE
#'
#' # To test if many things are metadata of an object
#' isMeta(c("age","groups"), myRNA) # FALSE, TRUE
#'
#' # 'return.values' input is especially useful in these cases.
#' isMeta(c("age","groups"), myRNA,
#'     return.values = TRUE)
#'
#' # Alternatively, to see all metadata slots of an object, use getMetas
#' getMetas(myRNA)
#'
#' @author Daniel Bunis
#' @export
#' @import ggplot2
#' @importFrom utils packageVersion

isMeta <- function(test, object, return.values=FALSE){
    if (return.values) {
        return(test[isMeta(test, object, return.values=FALSE)])
    } else {
        metas <- getMetas(object)
        if (is(object,"Seurat") || is(object,"seurat")) {
            metas <- c(metas,"ident")
        }
        return(test %in% metas)
    }
}

#### getMetas: prints the names of all the metadata lists for the object ####
#' Returns the names of all meta.data slots of a target object.
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param names.only Logical, \code{TRUE} by default, which sets whether just the names should be output versus the entire metadata dataframe.
#' @return A string vector of the names of all metadata slots of the \code{object}, or alternatively the entire dataframe of metadatas if \code{names.only} is set to \code{FALSE}
#' @seealso
#' \code{\link{isMeta}} for checking if certain metadata slots exist in an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#'
#' # To see all metadata slots of an object
#' getMetas(myRNA)
#'
#' # To retrieve the entire metadata matrix
#' getMetas(myRNA, names.only = FALSE)
#'
#' @author Daniel Bunis
#' @export

getMetas <- function(object, names.only = TRUE){
    metadata <-
        if (is(object,"SummarizedExperiment")) {
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
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param adjustment A recognized string indicating whether numeric metadata should be used directly (default) versus adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' 
#' Ignored if the target metadata is not numeric.
#' @param adj.fxn A function which takes a vector (of metadata values) and returns a vector of the same length.
#' 
#' For example, \code{function(x) \{log2(x)\}} or \code{as.factor}
#' @return A named vector.
#' @details
#' Retrieves the values of a metadata slot from \code{object}, or the clustering slot if \code{meta = "ident"} and the \code{object} is a Seurat.
#' 
#' If \code{adjustment} or \code{adj.fxn} are provided, then these requested adjustments are applied to these values (\code{adjustment} first).
#' Note: Alterations via \code{adjustment} are only applied when metadata is numeric, but \code{adj.fxn} alterations are applied to metadata of any type.
#' 
#' Lastly, outputs these values are named as the cells'/samples' names.
#' @seealso
#' \code{\link{metaLevels}} for returning just the unique discrete identities that exist within a metadata slot
#'
#' \code{\link{getMetas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{isMeta}} for testing whether something is the name of a metadata slot
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' meta("groups", object = myRNA)
#' 
#' myRNA$numbers <- seq_len(ncol(myRNA))
#' meta("numbers", myRNA, adjustment = "z-score")
#' meta("numbers", myRNA, adj.fxn = as.factor)
#' meta("numbers", myRNA, adj.fxn = function(x) \{log2(x)\})
#'
#' @author Daniel Bunis
#' @export

meta <- function(meta, object,
    adjustment = NULL, adj.fxn = NULL) {

    if (!isMeta(meta, object)) {
        stop(dQuote(meta)," is not a metadata of 'object'")
    }
    
    # Retrieve target metadata's values
    if (meta=="ident" && !is(object,"SummarizedExperiment")) {
        # Seurat clustering
        if (is(object, "Seurat")) {
            .error_if_no_Seurat()
            values <- Seurat::Idents(object)
        } else {
            values <- object@ident
        }
    } else {
        
        values <- getMetas(object, names.only = FALSE)[,meta, drop = TRUE]
    }

    # Add adjustments
    if (is.numeric(values)) {
        
        if (!is.null(adjustment) && !is.na(adjustment)) {
            if (adjustment=="z-score") {
                values <- as.numeric(scale(values))
            }
            if (adjustment=="relative.to.max") {
                values <- values/max(values)
            }
        }
    }
        
    if (!is.null(adj.fxn)) {
        values <- adj.fxn(values)
    }
    
    # Add names
    names(values) <- .all_cells(object)
    
    values
}

#' Gives the distinct values of a meta.data slot (or ident)
#'
#' @param meta quoted "meta.data.slot" name = REQUIRED. the meta.data slot whose potential values should be retrieved.
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param used.only TRUE by default, for target metadata that are factors, whether levels nonexistent in the target data should be ignored.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @return String vector, the distinct values of a metadata slot (factor or not) among all cells/samples, or for a subset of cells/samples.
#' (Alternatively, returns the distinct values of clustering if \code{meta = "ident"} and the object is a \code{Seurat} object).
#' @seealso
#' \code{\link{meta}} for returning an entire metadata slots of an \code{object}, not just the potential levels
#'
#' \code{\link{getMetas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{isMeta}} for testing whether something is the name of a metadata slot
#' @examples
#' example(importDittoBulk, echo = FALSE)
#'
#' metaLevels("clustering", object = myRNA)
#'
#' # Note: Set 'used.only' (default = TRUE) to FALSE to show unused levels
#' #  of metadata that are already factors.  By default, only the in use options
#' #  of a metadata are shown.
#' metaLevels("clustering", myRNA,
#'     used.only = FALSE)
#'
#' @author Daniel Bunis
#' @export

metaLevels <- function(meta, object, cells.use = NULL, used.only = TRUE){
    if (!isMeta(meta, object)) {
        stop(dQuote(meta)," is not a metadata of 'object'")
    }
    if (used.only) {
        meta.values <- as.character(meta(meta, object))
    } else {
        meta.values <- meta(meta, object)
    }
    if (!is.null(cells.use)) {
        all.cells <- .all_cells(object)
        cells.use <- .which_cells(cells.use, object)
        meta.values <- meta.values[all.cells %in% cells.use]
    }
    levels(as.factor(meta.values))
}
