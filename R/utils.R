#' @importFrom SingleCellExperiment reducedDim int_metadata
#' @importFrom SummarizedExperiment assay assays
#' @importFrom cowplot ggdraw get_legend

.remove_legend <- function(ggplot) {
    ggplot + theme(legend.position = "none")
}

.grab_legend <- function(ggplot) {
    cowplot::ggdraw(cowplot::get_legend(ggplot))
}

.all_cells <- function(object = DEFAULT) {
    if (is(object, "Seurat") || is(object,"SingleCellExperiment")) {
        # Seurat-v3, bulk or single-cell SCE
        return(colnames(object))
    } else {
        # Seurat-v2
        return(object@cell.names)
    }
}

.which_cells <- function(cells.use, object = DEFAULT) {
    all.cells <- .all_cells(object)
    if (is.null(cells.use)) {
        return(all.cells)
    }
    if (is.logical(cells.use)) {
        return(all.cells[cells.use])
    }
    # If here, not a logical or NULL, thus should already be a list of names
    cells.use
}

.var_OR_get_meta_or_gene <- function(var, object = DEFAULT,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL) {

    OUT <- var
    if (length(var)==1 && typeof(var)=="character") {
        if (isMeta(var, object)) {
            OUT <- meta(var, object)
        }
        if (isGene(var, object, assay)) {
            OUT <- gene(var, object, assay, slot, adjustment)
        }
    }
    names(OUT) <- .all_cells(object)
    OUT
}

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

.which_data <- function(
    assay = .default_assay(object), slot = .default_slot(object), object) {

    if (is(object,"SingleCellExperiment")) {
        return(SummarizedExperiment::assay(object, assay))
    }
    if (is(object,"Seurat")) {
        .error_if_no_Seurat()
        return(Seurat::GetAssayData(object, assay = assay, slot = slot))
    }
    if (is(object,"seurat")) {
        return(eval(expr = parse(text = paste0(
                "object@",slot))))
    }
}

.extract_Reduced_Dim <- function(reduction.use, dim=1, object) {
    # If object is a Seurat object
    if (is(object,"seurat")) {
        embeds <- eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@cell.embeddings")))
        colnames(embeds) <- paste0(eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@key"))),seq_len(ncol(embeds)))
    }
    if (is(object,"Seurat")) {
        .error_if_no_Seurat()
        embeds <- Seurat::Embeddings(object, reduction = reduction.use)
    }
    if (is(object,"SingleCellExperiment")) {
        embeds <- SingleCellExperiment::reducedDim(
            object, type = reduction.use, withDimnames=TRUE)
    }
    OUT <- list(embeds[,dim])
    OUT[2] <- colnames(embeds)[dim]
    names(OUT) <- c("embeddings","name")
    OUT
}

.gen_key <- function (reduction.use){
    key <- reduction.use
    if (grepl("pca|PCA", reduction.use)) {
        key <- "PC"
    }
    if (grepl("cca|CCA", reduction.use)) {
        key <- "CC"
    }
    if (grepl("cca.aligned", reduction.use)) {
        key <- "aligned.CC"
    }
    if (grepl("ica|ICA", reduction.use)) {
        key <- "IC"
    }
    if (grepl("tsne|tSNE|TSNE", reduction.use)) {
        key <- "tSNE_"
    }
    key
}

.make_hover_strings_from_vars <- function(
    data.hover, object, assay, slot, adjustment) {

    # Overall: if do.hover=TRUE and data.hover has a list of genes / metas called
    #   c(var1, var2, var3, ...), then for all cells, make a string:
    #   "var1: var1-value\nvar2: var2-value\nvar3: var3-value\n..."
    #   vars that are not genes or metadata are ignored.
    fillable <- vapply(
        seq_along(data.hover),
        function(i)
            (isMeta(data.hover[i],object) ||
                isGene(data.hover[i],object, assay)),
        logical(1))
    data.hover <- data.hover[fillable]
    if (is.null(data.hover)) {
        stop("No genes or metadata names added to `hover.data`")
    }

    # Create dataframe to contain the hover.info
    features.info <- data.frame(row.names = .all_cells(object))
    features.info <- vapply(
        data.hover,
        function(this.data)
            as.character(.var_OR_get_meta_or_gene(
                this.data,object, assay, slot, adjustment)),
        character(nrow(features.info)))
    names(features.info) <- data.hover[fillable]

    # Convert each row of dataframe to 'colname1: data1\ncolname2: data2\n...'
    hover.strings <- .make_hover_strings_from_df(features.info)
}

.make_hover_strings_from_df <- function(df){
    vapply(
        seq_len(nrow(df)),
        function(row){
            paste(as.character(vapply(
                seq_len(ncol(df)),
                function(col){
                    paste0(names(df)[col],": ",df[row,col])
                }, FUN.VALUE = character(1))
                ),collapse = "\n")
        }, FUN.VALUE = character(1))
}

.rename_and_or_reorder <- function(orig.data, reorder = NULL, relabels = NULL) {
    if (is.numeric(orig.data)) {
        return(orig.data)
    }
    rename.args <- list(x = orig.data)
    if (!(is.null(reorder))) {
        rename.args$levels <- levels(factor(rename.args$x))[reorder]
    }
    if (!(is.null(relabels))) {
        rename.args$labels <- relabels
    }
    do.call(factor, args = rename.args)
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

.is_bulk <- function(object) {
    OUT <- FALSE
    if (is(object,"SingleCellExperiment")) {
        if (!is.null(SingleCellExperiment::int_metadata(object)$bulk)) {
            if (SingleCellExperiment::int_metadata(object)$bulk) {
                OUT <- TRUE
            }
        }
    }
    OUT
}

.error_if_no_Seurat <- function() {
    if (!requireNamespace("Seurat")) {
        stop("Seurat installation required.")
    }
}

.error_if_no_plotly <- function() {
    if (!requireNamespace("plotly")) {
        stop("plotly installation required for using hover")
    }
}
