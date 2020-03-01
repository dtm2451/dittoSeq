#' import bulk sequencing data into a format that dittoSeq functions expect.
#' @rdname importDittoBulk
#' @param ... For the generic, additional arguments passed to specific methods.
#' @param x placeholder
#' @param reductions placeholder
#' @param metadata placeholder
#' @param combine_metadata placeholder
#' @return placeholder
#' One recommended assay to have, in addition to raw counts, is logcounts as this is what dittoSeq expression visualizations will look to by default.
#' A common way to generate one is through
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment "colData<-" colData rowData
#' @importFrom utils packageVersion
#' @importFrom SingleCellExperiment "int_metadata<-" int_metadata
#' @export
setGeneric("importDittoBulk", function(x, ...) standardGeneric("importDittoBulk"))

#' @rdname importDittoBulk
importDittoBulk.SummarizedExperiment <- function(
    x, reductions = NULL, metadata = NULL, combine_metadata = TRUE) {

    object <- as(x, "SingleCellExperiment")
    # Use SCE int_metadata to store dittoSeq version and that dataset is bulk
    int_metadata(object) <- c(int_metadata(object),
        dittoSeqVersion = packageVersion("dittoSeq"), bulk = TRUE)

    # Add metadata
    if (!is.null(metadata)) {
        if (combine_metadata) {
            # Add metadata from `x` to provided `metadata` dataframe.
            obj_metadata <- colData(object)
            dups <- colnames(obj_metadata) %in% colnames(metadata)
            # If any names are repeated, use from the provided `metadata`
            if (any(dups)) {
                message(
                    paste(colnames(obj_metadata)[dups], collapse = ", "),
                    " metadata originally within `x` was overwitten from provided `metadata`")
            }
            metadata <- cbind(obj_metadata[,!dups], metadata)
        }
        colData(object) <- metadata
    }

    # Add reductions
    if (!is.null(reductions)) {
        for (i in names(reductions)) {
            reducedDim(object, i) <- reductions[[i]]
        }
    }
    object
}

#' @rdname importDittoBulk
importDittoBulk.DGEList <- function(
    x, reductions = NULL, metadata = NULL, combine_metadata = TRUE) {

    ### Convert DGEList to Summarized Experiment while preserving as many
    ### optional slots as possible
    # Grab essential slots
    args <- list(assays = list(counts=x$counts),
                 colData = data.frame(samples = x$samples))
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
    se <- do.call(SummarizedExperiment, args)

    # Import as if the DGEList was always an SE
    importDittoBulk(se, reductions, metadata, combine_metadata)
}

