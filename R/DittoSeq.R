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

#### is.gene: Is this the name of a gene in my dataset? ####
#' Tests if input is the name of a gene in a target object.
#'
#' @param test String or vector of strings, the "potential.gene.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a gene within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed genes if \code{return.values = TRUE}.
#' @seealso
#' \code{\link{get.genes}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' get.genes(pbmc)
#'
#' # To test if something is a gene in an object:
#' is.gene("CD14", object = "pbmc") # TRUE
#' is.gene("CD12345", pbmc) # FALSE
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped.
#' DEFAULT <- "pbmc"
#' is.gene("CD14")
#'   # TRUE
#'
#' # To test if many things are genes of an object
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"))
#'
#' # return.values input is especially useful in these cases.
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"), return.values = TRUE)
#'
#' @export

is.gene <- function(test, object=DEFAULT, return.values = FALSE){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    if (return.values) {
        return(test[is.gene(test, object, return.values=FALSE)])
    } else {
        return(test %in% get.genes(object))
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

#### get.genes: prints the names of all the genes for a Seurat or RNAseq ####
#' Returns the names of all genes of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @return A string vector, returns the names of all genes of the \code{object}.
#' @seealso
#' \code{\link{is.gene}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' get.genes(pbmc)
#'
#' @export

get.genes <- function(object=DEFAULT){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    rownames(.which_data("normalized", object))
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


#### gene: for extracting the expression values of a particular gene for all cells/samples ####
#' Returns the values of a gene for all cells/samples
#'
#' @param gene               quoted "gene" name = REQUIRED. the gene whose expression data should be retrieved.
#' @param object             "name" of the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param data.type          Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), or "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation)? Default = "normalized"
#' @return Returns the values of a meta.data slot, or the ident (clustering) slot if "ident" was given and the object is a Seurat object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' gene("CD14", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' gene("CD14")
#' @export

gene <- function(gene, object=DEFAULT, data.type = "normalized"){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }

    # Recursive functions for non
    if (data.type == "relative") {
        return(as.numeric(scale(gene(gene, object, "normalized"))))
    }
    if (data.type == "normalized.to.max") {
        exp <- gene(gene, object, "normalized")
        return(exp/max(exp))
    }
    if (data.type == "raw.normalized.to.max") {
        exp <- gene(gene, object, "raw")
        return(exp/max(exp))
    }

    exp <- .which_data(data.type, object)[gene,]
    names(exp) <- .all_cells(object)

    exp
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
