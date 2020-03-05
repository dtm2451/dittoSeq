#' import bulk sequencing data into a format that dittoSeq functions expect.
#' @rdname importDittoBulk
#' @param ... For the generic, additional arguments passed to specific methods.
#' @param x a \code{\linkS4class{DGEList}}, or \code{\linkS4class{SummarizedExperiment}} (includes \code{DESeqDataSet}) class object containing the sequencing data to be imported
#' @param reductions a named list of dimensionality reduction embeddings matrices.
#' names will become the names of the dimensionality reductions and how each will be used with the \code{reduction.use} input of \code{dittoDimPlot}
#' rows of the matrices should represent the different cells/samples of the dataset, and columns the different dimensions
#' @param metadata a data.frame like object containing columns of extra information about the cells/samples (rows).
#' The names of these columns can then be used to tretrieve and plot such data in any dittoSeq visualizations.
#' @param combine_metadata Logical which sets whether original colData (DESeqDataSet/SummarizedExperiment) or $samples (DGEList) from x should be retained.
#' @return A \code{\linkS4class{SingleCellExperiment}} object containing all assays (DESeqDataSet or SummarizeedExperiment) or all common slots (DGEList) of the input \code{x},
#' as well as any dimensionality reductions provided to \code{reductions}, and any provided \code{metadata} stored in colData.
#'
#' When \code{combine_metadata} is set to FALSE, metadata inside \code{x} (colData or $samples) is ignored entirely.
#' When \code{combine_metadata} is TRUE (the default), metadata inside \code{x} is combined with what is provided to the \code{metadata} input; but names must be unique, so when there are similarly named slots, the values provided to the \code{metadata} input are used.
#'
#' @seealso \code{\linkS4class{SingleCellExperiment}} for more information about this storage system.
#'
#' @note One recommended assay to create if it is not already present in your dataset, is a log-normalized version of the counts data.
#' The logNormCounts function of the scater package is an easy way to make such a slot.
#' dittoSeq defaults to grabbing expression data from an assay named logcounts > normcounts > counts
#'
#' @examples
#' library(SingleCellExperiment)
#' nsamples <- 20
#' conditions <- factor(rep(1:4, each=5))
#' exp <- matrix(rpois(20000, 5), ncol=nsamples)
#' colnames(exp) <- paste0("sample", seq_len(ncol(exp)))
#' rownames(exp) <- paste0("gene", seq_len(nrow(exp)))
#' logexp <- log2(exp + 1)
#'
#' # Make this into a SummarizedExperiment
#' se <- SummarizedExperiment(
#'     list(counts = exp, logcounts = logexp),
#'     colData = data.frame(conditions))
#'
#' # import
#' myRNA <- importDittoBulk(se)
#'
#' ### For DESeq2, how we might have made this:
#' # library(DESeq2)
#' # dds <- DESeqDataSetFromMatrix(
#' #     exp, data.frame(conditions), ~ conditions)
#' # dds <- DESeq(dds)
#' # dds_ditto <- importDittoBulk(dds)
#'
#' ### For edgeR
#' # library(edgeR)
#' # dgelist <- DGEList(counts=exp, group=conditions)
#' # dge_ditto <- importDittoBulk(dgelist)
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom edgeR DGEList
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment "colData<-" colData rowData
#' @importFrom utils packageVersion
#' @importFrom SingleCellExperiment "int_metadata<-" int_metadata "reducedDim<-"
#' @importFrom S4Vectors DataFrame
#' @export
setGeneric("importDittoBulk", function(x, ...) standardGeneric("importDittoBulk"))

#' @rdname importDittoBulk
setMethod("importDittoBulk", "SummarizedExperiment", function(
    x, reductions = NULL, metadata = NULL, combine_metadata = TRUE) {

    object <- as(x, "SingleCellExperiment")
    # Use SCE int_metadata to store dittoSeq version and that dataset is bulk
    SingleCellExperiment::int_metadata(object) <- c(
        SingleCellExperiment::int_metadata(object),
        dittoSeqVersion = packageVersion("dittoSeq"),
        bulk = TRUE)

    # Add metadata
    if (!is.null(metadata)) {
        if (combine_metadata) {
            # Add metadata from `x` to provided `metadata` dataframe.
            obj_metadata <- SummarizedExperiment::colData(object)
            dups <- colnames(obj_metadata) %in% colnames(metadata)
            # If any names are repeated, use from the provided `metadata`
            if (any(dups)) {
                message(
                    paste(colnames(obj_metadata)[dups], collapse = ", "),
                    " metadata originally within 'x' was overwitten from provided 'metadata'")
            }
            metadata <- cbind(obj_metadata[,!dups, drop = FALSE], metadata)
        }
        SummarizedExperiment::colData(object) <- S4Vectors::DataFrame(metadata)
    }

    # Add reductions
    if (!is.null(reductions)) {
        if (is.null(names(reductions))) {
            stop("Elements of 'reductions' must be named to be added.")
        }
        for (i in names(reductions)) {
            if (i == "") stop("All elements of 'reductions' must be named.")
            SingleCellExperiment::reducedDim(object, i) <- reductions[[i]]
        }
    }
    object
})

#' @rdname importDittoBulk
setMethod("importDittoBulk", "DGEList", function(
    x, reductions = NULL, metadata = NULL, combine_metadata = TRUE) {

    ### Convert DGEList to Summarized Experiment while preserving as many
    ### optional slots as possible
    # Grab essential slots
    args <- list(
        assays = list(counts=x$counts),
        colData = data.frame(x$samples))
    # Add optional rowData
    rowData <- list()
    add_if_slot <- function(i, out = rowData) {
        if (!is.null(x[[i]])) {
            out <- c(out, list(x[[i]]))
            names(out)[length(out)] <- i
        }
        out
    }
    rowData <- add_if_slot("genes")
    rowData <- add_if_slot("AveLogCPM")
    rowData <- add_if_slot("common.dispersion")
    rowData <- add_if_slot("trended.dispersion")
    rowData <- add_if_slot("tagwise.dispersion")
    if (length(rowData)>0) {
        args$rowData <- data.frame(rowData)
    }
    # Add optional assay: offset
    if (!is.null(x$offset)) {
        args$assays <- list(counts = x$counts,
                            offset = x$offset)
    }
    # Make SE
    se <- do.call(SummarizedExperiment::SummarizedExperiment, args)

    # Import as if the DGEList was always an SE
    importDittoBulk(se, reductions, metadata, combine_metadata)
})

#' @rdname importDittoBulk
setMethod("importDittoBulk", "list", function(
    x, reductions = NULL, metadata = NULL) {

    # Check that the elements of x are named matrices with equal ncols
    ncol1 <- ncol(x[[1]])
    ncol.same <- all(vapply(
        seq_along(x),
        function (ind) {ncol(x[[ind]])==ncol1},
        FUN.VALUE = logical(1)))
    if (is.null(names(x)) || !ncol.same) {
        stop("Elements of 'x' should be named and should all have the same number of columns.")
    }

    # Create SummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(assays = x)

    # Import as if the DGEList was always an SE
    importDittoBulk(se, reductions, metadata, combine_metadata = TRUE)
})
