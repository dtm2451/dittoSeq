.all_cells <- function(object) {
    if (is(object, "Seurat") || is(object,"SingleCellExperiment")) {
        # Seurat-v3, bulk or single-cell SCE
        return(colnames(object))
    } else {
        # Seurat-v2
        return(object@cell.names)
    }
}

.which_cells <- function(cells.use, object) {
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

.var_OR_get_meta_or_gene <- function(var, object,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL) {

    OUT <- var
    cells <- .all_cells(object)
    if (length(var)==1 && is.character(var)) {
        if (isMeta(var, object)) {
            OUT <- meta(var, object)
        } else if (isGene(var, object, assay)) {
            OUT <- gene(var, object, assay, slot, adjustment)
        }
    }

    if (length(OUT)!=length(cells)) {
        stop("'var' is not a metadata or gene nor equal in length to ncol('object')")
    }
    names(OUT) <- cells
    OUT
}

.add_by_cell <- function(dat, target, name, object,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL, reorder = NULL, relabels = NULL, mult = FALSE) {

    # Extracts metadata or gene expression is target is length 1, then adds
    # this data to the dat dataframe as a column named name. If target length
    # = ncol(object), its values are used directly.
    # For discrete (factor) data, reorder and relabel are used to do just that.
    #
    # *In mult = TRUE mode, vector targets and vactor name are dealt with as
    # if there is 1 meta/gene and associated name per index of these vectors,
    # and reorder/delabels are not used.

    if (mult) {
        for (i in seq_along(target)) {
            dat <- .add_by_cell(dat, target[i], name[i], object, assay, slot,
                adjustment, reorder = NULL, relabels = NULL, mult = FALSE)
        }
    } else if (!is.null(target)) {
        # Obtain and reorder values
        values <- .var_OR_get_meta_or_gene(
            target, object, assay, slot, adjustment)
        values <- .rename_and_or_reorder(values, reorder, relabels)
        # Add
        dat <- cbind(dat, values)
        # Set the name
        names(dat)[ncol(dat)] <- name
    }

    dat
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
        if (is.null(colnames(embeds))) {
            colnames(embeds) <-
                paste0(.gen_key(reduction.use), seq_len(ncol(embeds)))
        }
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
