#### isGene: Is this the name of a gene in my dataset? ####
#' Tests if input is the name of a gene in a target object.
#'
#' @param test String or vector of strings, the "potential.gene.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a gene within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed genes if \code{return.values = TRUE}.
#' @seealso
#' \code{\link{getGenes}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' getGenes(pbmc)
#'
#' # To test if something is a gene in an object:
#' isGene("CD14", object = "pbmc") # TRUE
#' isGene("CD12345", pbmc) # FALSE
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped.
#' DEFAULT <- "pbmc"
#' isGene("CD14")
#'   # TRUE
#'
#' # To test if many things are genes of an object
#' isGene(c("CD14", "IL32", "CD3E", "CD12345"))
#'
#' # return.values input is especially useful in these cases.
#' isGene(c("CD14", "IL32", "CD3E", "CD12345"), return.values = TRUE)
#'
#' @author Daniel Bunis
#' @export

isGene <- function(test, object=DEFAULT, return.values = FALSE){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    if (return.values) {
        return(test[isGene(test, object, return.values=FALSE)])
    } else {
        return(test %in% getGenes(object))
    }
}

#### getGenes: prints the names of all the genes for a Seurat or RNAseq ####
#' Returns the names of all genes of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @return A string vector, returns the names of all genes of the \code{object}.
#' @seealso
#' \code{\link{isGene}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' getGenes(pbmc)
#'
#' @author Daniel Bunis
#' @export

getGenes <- function(object=DEFAULT){
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    rownames(.which_data(object=object))
}

#### gene: for extracting the expression values of a particular gene for all cells/samples ####
#' Returns the expression values of a gene for all cells/samples
#'
#' @param gene quoted "gene" name = REQUIRED. the gene whose expression data should be retrieved.
#' @param object "name" of the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param assay,slot single strings or integer that set which data to use.
#' Seurat and SingleCellExperiments deal with these differently, so be sure to check the documentation for whichever object you are using.
#' When not provided, these typical defaults for the provided \code{object} class are used:
#' \itemize{
#' \item{SingleCellExperiment (single-cell or bulk data): \code{assay} = "logcounts", "normcounts", "counts", or the first element of assays(object), \code{slot} not used}
#' \item{Seurat-v3: \code{assay} = DefaultAssay(object), \code{slot} = "data"}
#' \item{Seurat-v2: \code{assay} not used, \code{slot} = "data"}
#' }
#' @param adjustment Should expression data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @return Returns the expression values of a gene for all cells/samples.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' gene("CD14", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' gene("CD14")
#'
#' @author Daniel Bunis
#' @export

gene <- function(
    gene, object=DEFAULT,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = c(NULL,"z-score", "relative.to.max")){

    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    adjustment <- match.arg(adjustment)

    # Recursive functions for adjustments
    if (adjustment=="z-score") {
        return(as.numeric(scale(gene(gene,object,assay,slot))))
    }
    if (adjustment=="relative.to.max") {
        exp <- gene(gene,object,assay,slot)
        return(exp/max(exp))
    }

    exp <- as.vector(.which_data(assay, slot, object)[gene,])
    names(exp) <- .all_cells(object)

    exp
}

