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

