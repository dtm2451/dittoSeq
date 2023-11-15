#### isGene: Is this the name of a gene in my dataset? ####
#' Tests if input is the name of a gene in a target object.
#'
#' @param test String or vector of strings, the "potential.gene.name"(s) to check for.
#' @param assay single string or integer that sets which set of seq data inside the object to check.
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#'
#' @inheritParams gene
#'
#' @return Returns a logical vector indicating whether each instance in \code{test} is a rowname within the requested \code{assay} of the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed rownames if \code{return.values = TRUE}.
#' @seealso
#' \code{\link{getGenes}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#'
#' # To see the first 10 genes of an object of a particular assay
#' getGenes(myRNA, assay = "counts")[1:10]
#'
#' # To see all genes of an object for the default assay that dittoSeq would use
#' # leave out the assay input (again, remove `head()`)
#' head(getGenes(myRNA))
#'
#' # To test if something is a gene in an object:
#' isGene("gene1", object = myRNA) # TRUE
#' isGene("CD12345", myRNA) # FALSE
#'
#' # To test if many things are genes of an object
#' isGene(c("gene1", "gene2", "not-a-gene", "CD12345"), myRNA)
#'
#' # 'return.values' input is especially useful in these cases.
#' isGene(c("gene1", "gene2", "not-a-gene", "CD12345"), myRNA,
#'     return.values = TRUE)
#'
#' @author Daniel Bunis
#' @export

isGene <- function(
    test, object,
    assay = .default_assay(object),
    return.values = FALSE,
    swap.rownames = NULL) {

    object <- .swap_rownames(object, swap.rownames)
    
    if (return.values) {
        return(test[isGene(test, object, assay, return.values=FALSE)])
    } else {
        return(test %in% getGenes(object, assay))
    }
}

#' Returns the names of all genes of a target object.
#'
#' @param assay Single string or integer that sets which set of seq data inside the object to check.
#'
#' @inheritParams gene
#'
#' @return A string vector, returns the names of all genes of the \code{object} for the requested \code{assay}.
#' @seealso
#' \code{\link{isGene}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' getGenes(object = myRNA, assay = "counts")
#'
#' # To see all genes of an object for the default assay that dittoSeq would use
#' # leave out the assay input
#' getGenes(myRNA)
#'
#' # Seurat
#' # pbmc <- Seurat::pbmc_small
#' # # To see all genes of an object of a particular assay
#' # getGenes(pbmc, assay = "RNA")
#'
#' @author Daniel Bunis
#' @export

getGenes <- function(
    object, assay = .default_assay(object), swap.rownames = NULL) {
    
    object <- .swap_rownames(object, swap.rownames)
    
    if (length(assay)>1) {
        return(
            do.call(
                c,
                lapply(
                    seq_along(assay),
                    function(i) {
                        getGenes(object, assay[i])
                    })
            )
        )
    }
    
    rownames(.which_data(object=object,assay = assay))
}

#### gene: for extracting the expression values of a particular gene for all cells/samples ####
#' Returns the expression values of a gene for all cells/samples
#'
#' @param gene quoted "gene" name = REQUIRED. the gene whose expression data should be retrieved.
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param assay,slot single strings or integers (SCEs and SEs) or a vector of such values that set which expression data to use.
#' See \code{\link{GeneTargeting}} for specifics and examples -- Seurat and SingleCellExperiment objects deal with these differently, and functionality additions in dittoSeq have led to some minimal divergence from the native methodologies.
#' @param adjustment Should expression data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param adj.fxn A function which takes a vector (of metadata values) and returns a vector of the same length.
#' 
#' For example, \code{function(x) \{log2(x)\}} or \code{as.factor}
#' @param swap.rownames optionally named string or string vector.
#' For SummarizedExperiment or SingleCellExperiment objects, its value(s) specifies the column name of rowData(object) to be used to identify features instead of rownames(object).
#' When targeting multiple modalities (alternative experiments), names can be used to specify which level / alternative experiment (use 'main' for the top-level) individual values should be used for.
#' See \code{\link{GeneTargeting}} for more specifics and examples.
#' @return Returns the expression values of a gene for all cells/samples.
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' gene("gene1", object = myRNA, assay = "counts")
#'
#' # z-scored
#' gene("gene1", object = myRNA, assay = "counts",
#'     adjustment = "z-score")
#' 
#' # Log2'd
#' gene("gene1", object = myRNA, assay = "counts",
#'     adj.fxn = function(x) \{log2(x)\})
#'
#' # To see expression of the gene for the default assay that dittoSeq would use
#' # leave out the assay input
#' # (For this object, the default assay is the logcounts assay)
#' gene("gene1", myRNA)
#'
#' # Seurat (raw counts)
#' if (!requireNamespace("Seurat")) {
#'     gene("CD14", object = Seurat::pbmc, assay = "RNA", slot = "counts")
#' }
#'
#' @author Daniel Bunis
#' @export

gene <- function(
    gene, object,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL, adj.fxn = NULL,
    swap.rownames = NULL) {
    
    object <- .swap_rownames(object, swap.rownames)

    if (!isGene(gene, object, assay)) {
        stop(dQuote(gene)," is not a gene/feature of the targeted assay(s) of 'object'")
    }
    
    # Retrieve target values
    exp <- as.vector(.which_data(assay, slot, object)[gene,])
    
    # Add adjustments
    if (!is.null(adjustment) && !is.na(adjustment)) {
        if (adjustment=="z-score") {
            exp <- as.numeric(scale(exp))
        }
        if (adjustment=="relative.to.max") {
            exp <- exp/max(exp)
        }
    }
        
    if (!is.null(adj.fxn)) {
        exp <- adj.fxn(exp)
    }

    # Add names
    names(exp) <- .all_cells(object)

    exp
}
