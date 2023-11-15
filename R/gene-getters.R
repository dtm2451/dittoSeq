#### isGene: Is this the name of a gene in my dataset? ####
#' Tests if input is the name of a gene in a target object.
#'
#' @inheritParams gene
#' @param test String or vector of strings, the "potential.gene.name"(s) to check for.
#' @param assay single string or integer that sets which set of seq data inside the object to check.
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
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
#' @inheritParams gene
#' @param assay Single string or integer that sets which set of seq data inside the object to check.
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
#' 
#' @param adj.fxn A function which takes a vector (of metadata values) and returns a vector of the same length.
#' 
#' For example, \code{function(x) \{log2(x)\}} or \code{as.factor}
#' @param swap.rownames String. For SummarizeedExperiment or SingleCellExperiment objects, the column name of rowData(object) to be used to identify features instead of rownames(object).
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

#' @name GeneTargeting
#' @title Control of Gene/Feature targeting
#' @section Overview:
#' As of dittoSeq version 1.15.2, we made it possible to target genes / features from across multiple modalities.
#' Here, we describe intricacies of how 'assay', 'slot', and 'swap.rownames' inputs now work to allow for this purpose.
#'
#' Control of gene/feature targeting in dittoSeq functions aims to blend seamlessly with how similar control works in Seurat, SingleCellExperiment (SCE), and other packages that deal with these data structures.
#' However, as we've built in new features into dittoSeq, and the Seurat and SCE-package maintainers extend their tools as well, some divergence was to be expected.
#'
#' The way Seurat and SingleCellExperiment objects hold data from multiple modalities is quite distinct, thus it is worth describing each distinctly.
#'
#' It's also important to note, that \emph{both structures utilize the term 'assay', but they utilize it for distinct meanings.}
#' Keep that in mind because we chose to stick with the native terminologies within dittoSeq in order maintain intuitiveness with other Seurat or SCE data accession methods.
#' In other words, rather than enforcing a new consistent paradigm, the native Seurat 'assay' meaning is respected for Seurat objects, and the native SCE 'assay' meaning is respected for SCE objects.
#'
#' @section Control of Gene/Feature targeting in Seurat Objects:
#' For Seurat objects, dittoSeq uses of its \code{assay} and \code{slot} inputs for gene/feature retrieval control, and ultimately makes use of Seurat's GetAssayData function for extracting data. (See: '?SeuratObject::GetAssayData')
#'
#' To allow targeting of features across multiple modalities, we allow provision of multiple assay names to dittoSeq's version of the 'assay' input.
#' Internally, dittoSeq will then loop through all values of 'assay', making a separate calls to GetAssayData for each assay.
#'
#' Otherwise, dittoSeq's \code{assay} and \code{slot} inputs work exactly the same as described in Seurat's documentation.
#'
#' Phrased another way, it works via inputs:
#' \itemize{
#' \item \code{assay} - takes the name(s) of Seurat Assays to target. Examples: \code{"RNA"} or \code{c("RNA", "ADT")}
#' \item \code{slot} - "counts", "data", or "scale.data". Directs which 'slot' of data from the targeted assays to extract from. Example: \code{"data"}
#' }
#'
#' As an example, if you wanted to plot raw counts data from 1) the CD4 gene of the RNA assay and 2) the CD4.1 marker of an ADT assay, you would:
#' \itemize{
#' \item 1. point the \code{var} or \code{vars} input of the plotter to \code{c("CD4", "CD4.1")}
#' \item 2. target both modalities via \code{assay = c("RNA", "ADT")} (Note that "RNA" and "ADT" are the default assay names typically used, but you do need to match with what is in your own Seurat object if your assays are named differently.)
#' \item 3. target the raw counts data via \code{slot = "counts"}
#' }
#'
#' @section Control of Gene/Feature targeting in SingleCellExperiment Objects:
#' For SCE objects, dittoSeq makes use of its \code{assay} input for both modality and data form (the meaning of 'assay' for SCEs) control, and ultimately makes use of the \code{\link[SummarizedExperiment]{assay}} and \code{\link[SingleCellExperiment]{altExp}} functions for extracting data.
#'
#' Additionally, we allow use of the \code{swap.rownames} input to allow targeting & display of alternative gene/feature names. The implementation here is that rownames of the extracted assay data are swapped out for the given \code{\link[SummarizedExperiment]{rowData}} column of the object (or altExp).
#' When used, note that you will need to use these swapped names for targeting genes / features with \code{gene}, \code{var}, or \code{vars} inputs.
#'
#' \strong{In SCE objects} themselves, the primary modality's expression data are stored in 'assay's of the SCE object.
#' You might have one assay containing raw data, and another containing log-normalized data.
#' Additional details of genes/features of this modality, possibly including alternative gene names, can be stored in the object's 'rowData' slot.
#' When additional modalities are collected, the way to store them is via a nested SCE object called an "alternative experiment".
#' Any number of these can be stored in the 'altExps' slot of the SCE object.
#' Each alternative experiment can contain any number of assays.
#' Again each will often have one representing raw data and another representing a normalized form of that data.
#' And, these alternative experiments might also make use of their rowData to store additional characteristics or names of each gene/feature.
#'
#' \emph{The system feels a bit more complicated here, because the SCE system is itself a bit more complicated. But the hope is that this system becomes simple to work with once learned!}
#'
#' To allow targeting of features across multiple modalities, dittoSeq's \code{assay} input can be given:
#' \itemize{
#' \item Simplest form: a single string or string vector where values are either the names of an assay of the primary modality OR the name of an alternative experiment to target, with 'main' as an indicator for the primary modality and 'altexp' as a shortcut for indicating "the first altExp".
#' In this form, when 'main', 'altexp', or the actual name of an alternative experiment are used, the first assay of that targeted modality will be used.
#' \item Explicit form: a named string or named vector of string values where names indicate the modality/experiment to target and values indicate what assays of those experiments to target. Here again, you can use 'main' or 'altexp' as names to mean the primary modality and "the first altExp", respectively.
#' \item These methods can also be combined. A few examples:
#' \itemize{
#' \item Using the simplified method only: \code{assay = c('main', 'altexp', 'hto')} will target the first assays each of the main object, of the first alternative experiment of the object, and also of an alternative experiment named 'hto'.
#' \item Using the explicit form only: \code{assay = c('main'='logexp', 'adt'='clr', 'altexp'='raw')} will target 1) the logexp-named assay of the main object, 2) the clr-named assay of an alternative experiment named 'adt', and 3) the raw-named assay of the first alternative experiment of the object.
#' \item Using a combination of the two: \code{assay = c('logexp', 'adt'='clr')} will target 1) the logexp-named assay of the primary modilty, unless there is an alternative experiment named 'logexp' which will lead to grabbing the first assay of that modality, and 2) the clr-named assay of an alternative experiment named 'adt'.
#' }
#' }
#'
#' The \code{swap.rownames} input allows swapping to alternative names for genes/features via provision of a column name of rowData(object).
#' The values of that rowData column are then used to identify and label features of the moadilty's assays instead of the original rownames of the assays.
#' To allow swap.rownames to also work with the multi-modality access system in the most simplified way, the swap.rownames input also has both a simple and an explicit provision system:
#' \itemize{
#' \item Simple form: a single string or string vector where any accessed modality will be checked for the presence of a these values in colnames of their rowData. Whenever found, those columns will be used.
#' \item Explicit form: similar to the explicit assay use, a named string or named vector of string values where names indicate the modality/experiment to target and values indicate columns to look for among the given modality's rowData.
#' \item Examples:
#' \itemize{
#' \item Simplified1: Using \code{assay = c('main', 'altexp'), swap.rownames = "SYMBOL"} with an object where the primary modality rowData has a SYMBOL column and the first alternative experiment's rowData is empty, will lead to swapping to the SYMBOL values for main modality features and use of original rownames for the alternative experiment's features.
#' \item Simplified2: Using \code{assay = c('main', 'altexp'), swap.rownames = "SYMBOL"} with an object where both modalities' rowData have a SYMBOL column, will lead to swapping to the SYMBOL values both modalities.
#' \item Explicit: Using \code{assay = c('main', 'altexp'), swap.rownames = c(main="SYMBOL")} with an object where both modalities' rowData have a SYMBOL column, will lead to swapping to the SYMBOL values for main modality only.
#' }
#' }
#'
#' As a full example, if you wanted to plot from 1) the raw 'counts' assay for a CD4 gene of the primary modality and 2) the normalized 'logexp' assay for a CD4.1 marker of an alternative experiment assay named 'ADT', but where 3) the rownames of these modalities are Ensembl ids while gene symbol names are held in a rowData column of both modalities that is named "symbols", the simplest provision method is:
#' \itemize{
#' \item 1. point the \code{var} or \code{vars} input of the plotter to \code{c("CD4", "CD4.1")}
#' \item 2. target the counts assay of the primary modality and logexp assay of the ADT alternative experiment via \code{assay = c('counts', ADT = 'logexp')}
#' \item 3. swap to the symbol names of features from both modlities by also giving \code{swap.rownames = "symbols"}
#' }
#' @author Dan Bunis
# NULL