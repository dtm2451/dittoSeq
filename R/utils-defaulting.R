.preferred_or_first <- function(opts, prefs) {
    # Given string vectors opts and prefs,
    # with prefs being a set of targets in decreasing preferrence order and all
    # lower case,
    # outputs the element of opts that contains the earliest element of prefs.
    # If no options contain an element of prefs, outputs the first opts.

    out <- opts[1]
    # tolower is used to allow upper/lower case to be ignored
    preferred <- match(prefs, tolower(opts))
    if (any(!is.na(preferred))) {
        out <- (opts[preferred[!is.na(preferred)]])[1]
    }
    out
}

.leave_default_or_null <- function(
    target, default, null.if = FALSE, default.when = "make") {
    # Handles much of dittoSeq's titles defaulting process
    # Takes in 'target' and outputs:
    #  - 'default' string when 'target' == 'default.when'
    #  - NULL when logical provided to 'null.if' is TRUE.
    #  - 'target' otherwise
    if (!is.null(target)) {
        if (identical(target,default.when)) {
            if (null.if) {
                target <- NULL
            } else {
                target <- default
            }
        }
    }
    target
}

#' @importFrom SummarizedExperiment assays
.default_assay <- function(object) {
    # Decides which assay should be default on prefs or defaults
    if (is(object, "SummarizedExperiment")) {
        # prefer logcounts > normcounts > counts > first assay
        return(.preferred_or_first(
            names(SummarizedExperiment::assays(object)),
            c("logcounts","normcounts","counts")))
    }
    if (is(object, "seurat")) {
        # no assays for Seurat-v2
        return(NA)
    }
    if (is(object, "Seurat")) {
        # use default assay
        .error_if_no_Seurat()
        return(Seurat::DefaultAssay(object))
    }
}

.default_assay_raw <- function(object) {
    # Decides which assay should be default on prefs or defaults
    if (is(object, "SummarizedExperiment")) {
        # prefer logcounts > normcounts > counts > first assay
        return(.preferred_or_first(
            names(SummarizedExperiment::assays(object)),
            c("counts","logcounts","normcounts")))
    }
    if (is(object, "seurat")) {
        # no assays for Seurat-v2
        return(NA)
    }
    if (is(object, "Seurat")) {
        # use default assay
        .error_if_no_Seurat()
        return(Seurat::DefaultAssay(object))
    }
}

.default_slot <- function(object) {
    # Decides what slot should be by default
    if (is(object, "SummarizedExperiment")) {
        # no slots for SCEs
        return(NA)
    } else {
        # default to the normalized data for Seurats
        return("data")
    }
}

.default_slot_raw <- function(object) {
    # Decides what slot should be by default
    if (is(object, "SummarizedExperiment")) {
        # no slots for SCEs
        return(NA)
    } else {
        # default to the normalized data for Seurats
        return("counts")
    }
}

.default_order <- function(object, annot.by) {
    # Sets the default for dittoHeatmap's 'order.by' to either NULL (no
    # ordering), or the first element of 'annot.by' if the object is an SCE and
    # annot.by has values.
    if (!is.null(annot.by) && !isBulk(object)) {
        return(annot.by[1])
    } else {
        return(NULL)
    }
}

#' @importFrom SingleCellExperiment reducedDim
.default_reduction <- function(object) {
    # Sets the default for dittoDimPlot's 'reduction.use' input
    # Capitalization-ignored, prefers umap > tsne > pca > the first reduciton.
    opts <- getReductions(object)
    if (is.null(opts)) {
        stop("No dimensionality reduction slots in 'object'. Add one, or provide embeddings directly to 'reduction.use'.")
    }
    use <- .preferred_or_first(opts, c("umap","tsne","pca"))
    use
}
