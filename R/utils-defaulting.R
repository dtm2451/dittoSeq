.preferred_or_first <- function(opts, prefs) {
    # Given a String vectors options and prefs,
    # with prefs being a set of targets in decreasing preferrence order and all
    # lower case, outputs the first element of options that is found to be an
    # element of prefs.
    # tolower is used to allow upper/lower case to be ignored
    # If no options contain an element of prefs, outputs the first option.
    out <- opts[1]
    preferred <- match(prefs, tolower(opts))
    if (any(!is.na(preferred))) {
        out <- (opts[preferred[!is.na(preferred)]])[1]
    }
    out
}

.leave_default_or_null <- function(
    target, default, null.if = FALSE, default.when = "make") {
    if (!is.null(target)) {
        if (target==default.when) {
            if (null.if) {
                target <- NULL
            } else {
                target <- default
            }
        }
    }
    target
}

.default_assay <- function(object) {
    #Decide which assy should be used
    if (is(object, "SingleCellExperiment")) {
        return(.preferred_or_first(
            names(SummarizedExperiment::assays(object)),
            c("logcounts","normcounts","counts")))
    }
    if (is(object, "seurat")) {
        return(NA)
    }
    if (is(object, "Seurat")) {
        .error_if_no_Seurat()
        return(Seurat::DefaultAssay(object))
    }
}

.default_slot <- function(object) {
    #Decide which assy should be used
    if (is(object, "SingleCellExperiment")) {
        return(NA)
    } else {
        return("data")
    }
}

.default_order <- function(object, annotation.metas) {
    if (!is.null(annotation.metas) && !.is_bulk(object)) {
        return(annotation.metas[1])
    } else {
        return(NULL)
    }
}

.default_reduction <- function(object) {
    # Use umap > tsne > pca, or whatever the first reduction slot is.
    opts <- getReductions(object)
    if (is.null(opts)) {
        stop("No dimensionality reduction slots in 'object'")
    }
    use <- .preferred_or_first(opts, c("umap","tsne","pca"))
    use
}
