#' Creates an RNAseq object from a DESeq2 processed bulk RNAseq data object.
#'
#' @description The first step of visualization of DESeq-analyzed bulk RNAseq data with dittoSeq is running this function to transform the data into the \linkS4class{RNAseq} structure that dittoSeq functions expect.
#' @param dds The output of running DESeq() on your data. = The DESeq2 object for your data. REQUIRED.
#' @param normalized.data (optional) a matrix containing normalized data for genes (rows) and samples (columns) in the dds.
#' column names must be the sample names, and row names must be gene names.
#' When provided, this matrix is used directly and not re-calculated.
#' @param blind Logical which sets whether rlog estimation should be blinded to sample info. Run \code{\link[DESeq2]{rlog}} for more info about whether it should.
#' @param counts (optional) a matrix containing raw counts data.
#' When provided, this matrix is used for the RNAseq \code{counts} slot. Otherwise, counts matrix embedded within the \code{dds} is used.
#' @param reductions (optional) for advanced users, a named list containing individual dimensionality reductions perviously run.
#' @return Outputs an \code{\linkS4class{RNAseq}} object.
#' @note dittoSeq is built to be a companion package for aiding analysis of RNAseq datasets of diverse structures.
#' It works quite well for this purpose, but is not meant to be used as the main tool for data-processing and normalization.
#' Inclusion of automatic normalization of bulk RNAseq data in import functions is provided for ease-of-use, but it is the responsibility of the user to provide properly normalized data
#' (to the `normalized.data` inputs of import functions) when the default settings for rlog (DESeq2) or cpm (edgeR/Limma-Voom) calculation are improper.
#' @details This function takes in a \code{dds} DESeqDataSet (the output of running \code{DESeq()})
#' and generates an \code{\linkS4class{RNAseq}} object with
#' raw \code{counts}, normalized \code{data}, \code{samples} names, and \code{meta.data} slots automatically pulled and filled in from the \code{dds}.
#' This effectively converts the data into the structure expected by all dittoSeq visualization and helper functions.
#'
#' For the \code{data} slot, the user must either provide a pre-calculated normalized data matrix or determine whether automatic DESeq2 regularized log calculation should be blinded.
#' For the pre-calculated data option, user should supply this directly as a matrix to the \code{normalized.data} input.
#' For the automatic calculation of regularized log normalization option, supply either \code{TRUE} or \code{FALSE} to the \code{blind} input which will be passed through to the DESeq2 \code{\link[DESeq2]{rlog}} function.
#' For more information on that, see \code{\link[DESeq2]{rlog}}.
#'
#' Metadata are pulled from the supplied \code{dds}. Then "Sample" names and "Nreads" metadata are generated automatically if not already included.
#' @seealso
#' \code{\link[DESeq2]{DESeq}} and \code{\link[DESeq2]{rlog}} for information about DESeq and regularized log calculation.
#'
#' \code{\link{addMetaRNAseq}} for how to add additional metadata to an \code{RNAseq} object
#'
#' \code{\link{addPrcomp}} and \code{\link{addDimReduction}} for how to add dimensionality reductions to an \code{\linkS4class{RNAseq}} object for use in
#' \code{\link{dittoDimPlot}} visualizations.
#' @examples
#'
#' # Generate mock RNAseq counts and a DESeq object from the mock data
#' # count tables from RNA-Seq data
#' counts.table <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#' colnames(counts.table) <- paste0("Sample",1:10)
#' rownames(counts.table) <- paste0("Gene",1:100)
#' conditions <- factor(rep(1:2, each=5))
#'
#' # object construction
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(
#'     counts.table, DataFrame(conditions), ~ conditions)
#' dds <- DESeq(dds)
#'
#' # Import
#' myRNA <- importDESeq2(dds, blind = FALSE)
#' @author Daniel Bunis
#' @export

importDESeq2 <- function(
    dds, normalized.data = NULL, blind, counts = NULL, reductions = NULL) {

    object <- new(
        "RNAseq", DEobject = dds, DEtype = "DESeq2")

    if (is.null(counts)) {
        object@counts <- DESeq2::counts(dds)
    } else {
        object@counts <- counts
    }

    if (!is.null(reductions)) {
        object@reductions <- reductions
    }

    object@samples <- colnames(object@counts)

    object@meta.data <- as.data.frame(SummarizedExperiment::colData(dds))
    if (is.null(object@meta.data$Sample)) {
        object@meta.data$Sample <- object@samples
    }
    if (is.null(object@meta.data$Nreads)) {
        object@meta.data$Nreads <- colSums(object@counts)
    }

    if (is.null(normalized.data)) {
        object@data <- SummarizedExperiment::assay(
            DESeq2::rlog(dds, blind = blind))
    } else {
        # Check that row and column names include all samples and gene names
        #   in the counts matrix, then add.
        if (!all(object@samples %in% colnames(normalized.data))) {
            stop("normalized.data must contain data for all samples")
        }
        if (!all(rownames(object@counts) %in% rownames(normalized.data))) {
            warning("Gene names in counts are missing in the provided normalized.data")
        }
        object@data <- normalized.data
    }

    object
}

#' Creates an RNAseq object from edgeR or limma-voom processed bulk RNAseq data objects.
#'
#' @description The first step of visualization of edgeR- or Limma-Voom-analyzed bulk RNAseq data with dittoSeq is running this function to transform the data into the \linkS4class{RNAseq} structure that dittoSeq functions expect.
#' @param DGEList a \link[edgeR]{DGEList-class} object, the main edgeR / Limma-Voom object for storing your data. REQUIRED.
#' @param normalization.fxn String which sets the counts data normalization function to use when calculating the \code{data} slot.
#' This is ignored when a matrix is provided to the \code{normalized.data} input.
#' @param normalized.data (optional) a matrix containing normalized data for all genes (rows) and samples (columns) in the DGEList.
#' When provided, this matrix is used directly and the \code{normalization.fxn} is not used to calculate the \code{data} slot from the \code{DGEList$counts} data.
#' @param counts (optional) a matrix containing raw counts data for all genes (rows) and samples (columns) in the DGEList.
#' When provided, this matrix is used, instead of the counts matrix embedded within the \code{DGEList}, for the RNAseq \code{counts} slot, but NOT for automatic calculation of the \code{data} slot.
#' @param reductions (optional) for advanced users, a named list containing individual dimensionality reductions perviously run.
#' @return Outputs an \code{\linkS4class{RNAseq}} object.
#' @note dittoSeq is built to be a companion package for aiding analysis of RNAseq datasets of diverse structures.
#' It works quite well for this purpose, but is not meant to be used as the main tool for data-processing and normalization.
#' Inclusion of automatic normalization of bulk RNAseq data in import functions is provided for ease-of-use, but it is the responsibility of the user to provide properly normalized data
#' (to the `normalized.data` inputs of import functions) when the default settings for rlog (DESeq2) or cpm (edgeR/Limma-Voom) calculation are improper.
#' @details This function takes in a \code{DGEList} data structure and generates a dittoSeq \code{\linkS4class{RNAseq}} object with
#' raw \code{counts}, normalized \code{data}, \code{samples} names, and \code{meta.data} slots automatically pulled and filled in from the \code{DGEList}.
#' This effectively converts the data into the structure expected by all dittoSeq visualization and helper functions.
#'
#' Because a common storage system for normalized counts within the DGEList structure is not set, for the \code{data} slot, this data is calculated with the edgeR \code{\link[edgeR]{cpm}} function.
#' Alternatively, this data can be directly supplied as a matrix to the \code{normalized.data} input, or
#' user's can provide the name of an alternative data normalization function to the \code{normalization.fxn} input to have the slot generated with this function instead (with default settings).
#'
#' Metadata are pulled from the supplied \code{DGEList$samples} slot. Then "Sample" names and "Nreads" metadata are generated automatically if not already included.
#' @seealso
#' \code{\link[edgeR]{edgeR}}, \code{\link[edgeR]{DGEList-class}}, and \code{\link[edgeR]{cpm}} for information about edgeR, the DGEList class, and the log2 counts per million calculation.
#'
#' \code{\link{addMetaRNAseq}} for how to add additional metadata to an \code{RNAseq} object
#'
#' \code{\link{addPrcomp}} and \code{\link{addDimReduction}} for how to add dimensionality reductions to an \code{\linkS4class{RNAseq}} object for use in
#' \code{\link{dittoDimPlot}} visualizations.
#' @examples
#'
#' # Generate mock RNAseq counts and a DGEList object from the mock data
#' # count tables from RNA-Seq data
#' counts.table <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#' colnames(counts.table) <- paste0("Sample",1:10)
#' rownames(counts.table) <- paste0("Gene",1:100)
#' conditions <- factor(rep(1:2, each=5))
#' # object construction
#' library(edgeR)
#' DGE <- DGEList(counts=counts.table,group=conditions)
#' DGE <- calcNormFactors(DGE)
#' design <- model.matrix(~conditions)
#' DGE <- estimateDisp(DGE,design)
#'
#' # Import
#' myRNA <- importEdgeR(DGE)
#' @author Daniel Bunis
#' @export

importEdgeR <- function(
    DGEList, normalized.data = NULL, normalization.fxn = "edgeR::cpm",
    counts = NULL, reductions = NULL) {

    object <- new(
        "RNAseq", DEobject = DGEList, DEtype = "EdgeR")

    if (is.null(counts)) {
        object@counts <- DGEList$counts
    } else {
        object@counts <- counts
    }

    if (!is.null(reductions)) {
        object@reductions <- reductions
    }

    object@samples <- colnames(object@counts)

    object@meta.data <- as.data.frame(DGEList$samples)
    if (is.null(object@meta.data$Sample)) {
        object@meta.data$Sample <- object@samples
    }
    if (is.null(object@meta.data$Nreads)) {
        object@meta.data$Nreads <- colSums(object@counts)
    }

    if (is.null(normalized.data)) {
        object@data <- eval(expr = parse(text =
            paste0(normalization.fxn,"(DGEList)")))
    } else {
        # Check that row and column names include all samples and gene names
        #   in the counts matrix, then add.
        if (!all(object@samples %in% colnames(normalized.data))) {
            stop("normalized.data must contain data for all samples")
        }
        if (!all(rownames(object@counts) %in% rownames(normalized.data))) {
            warning("Gene names in counts are missing in the provided normalized.data")
        }
        object@data <- normalized.data
    }

    object
}
