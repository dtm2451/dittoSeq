.all_cells <- function(object) {
    # Retrieves all cell/sample names / barcodes from a dataset
    if (is(object, "Seurat") || is(object,"SummarizedExperiment")) {
        # Seurat-v3, bulk or single-cell SCE
        return(colnames(object))
    } else {
        # Seurat-v2
        return(object@cell.names)
    }
}

.which_cells <- function(cells.use, object) {
    # converts a 'cells.use' given as string, logical, or numeric vector
    # into the string vector format expected internally by dittoSeq functions
    all.cells <- .all_cells(object)
    
    if (is.null(all.cells)) {
        stop(
            "Given 'object' has no cell/column names, which will inhibit dittoSeq's data gathering functions. ",
            "For most compatible 'object's, you can correct the issue with:\n    ",
            "`colnames(<object>) <- paste0('cell', seq_len(ncol(<object>)))`")
    }

    if (is.null(cells.use)) {
        # Returns all cells when 'cells.use' is NULL
        return(all.cells)
    }

    if (is.logical(cells.use)) {
        if (length(cells.use)!=length(all.cells)) {
            stop("'cells.use' length must equal the number of cells/samples in 'object' when given in logical form")
        }
        return(all.cells[cells.use])
    }

    if (is.numeric(cells.use)) {
        return(all.cells[cells.use])
    }

    cells.use
}

#' @importFrom SummarizedExperiment assay
.which_data <- function(
    assay = .default_assay(object), slot = .default_slot(object), object) {
    # Retrieves the required counts data from 'object'
    
    if (length(assay)>1) {
        return(
            do.call(
                rbind,
                lapply(
                    seq_along(assay),
                    function(i) {
                        .which_data(assay[i], slot, object)
                    })
            )
        )
    }

    if (is(object,"SummarizedExperiment")) {
        if (is(object,"SingleCellExperiment")) {
            ### altExp compatibility
            # Note: name='' check because that's given to unnamed elements of a half-named c() call.
            if (identical(NULL, names(assay)) || identical('', names(assay))) {
                # Simplified method: check 'altexp' & 'main' tokens, then top-level assay names, then altExp names
                if (assay=='altexp') {
                    return(SummarizedExperiment::assay(SingleCellExperiment::altExp(object)))
                } else if (assay=='main') {
                    return(SummarizedExperiment::assay(object))
                } else if (assay %in% SummarizedExperiment::assayNames(object)) {
                    return(SummarizedExperiment::assay(object, assay))
                } else {
                    return(SummarizedExperiment::assay(SingleCellExperiment::altExp(object, assay)))
                }
            } else {
                # Explicit (named) method: check 'altexp' & 'main' tokens, then explicit altExp names
                if (names(assay)=='altexp') {
                    return(SummarizedExperiment::assay(SingleCellExperiment::altExp(object), assay))
                } else if (names(assay)=='main') {
                    return(SummarizedExperiment::assay(object, assay))
                } else {
                    return(SummarizedExperiment::assay(SingleCellExperiment::altExp(object, names(assay)), assay))
                }
            }
        }
        return(SummarizedExperiment::assay(object, assay))
    }
    if (is(object,"Seurat")) {
        .error_if_no_Seurat()
        if (packageVersion("Seurat")>='5.0') {
            return(object[[assay]][slot])
        }
        return(Seurat::GetAssayData(object, assay = assay, slot = slot))
    }
    if (is(object,"seurat")) {
        return(eval(expr = parse(text = paste0(
                "object@",slot))))
    }
}

.var_OR_get_meta_or_gene <- function(var, object,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL) {
    # Turns 'var' strings refering to genes or metadata into their associated data
    # Otherwise, returns 'var' with cellname names added.

    OUT <- var
    
    cells <- .all_cells(object)
    object <- .swap_rownames(object, swap.rownames)
    
    if (length(var)==1 && is.character(var)) {
        if (isMeta(var, object)) {
            OUT <- meta(var, object)
        } else if (isGene(var, object, assay)) {
            OUT <- gene(var, object, assay, slot, adjustment)
        }
    }

    if (length(OUT)!=length(cells)) {
        stop(
            ifelse(length(var)==1, var, 'var'),
            " is not a gene of the targeted assay(s), a metadata, nor equal in length to ncol('object')")
    }
    names(OUT) <- cells
    OUT
}

.add_by_cell <- function(df = NULL, target, name, object,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL, reorder = NULL, relabels = NULL,
    mult = FALSE) {

    # Extracts metadata or gene expression if 'target' is the name of one,
    # or if length('target') = ncol(object), its values are used directly.
    # These values are added to the 'df' dataframe as a column named 'name'.
    #
    # For gene data, 'assay', 'slot', and 'adjustment' control how the data is
    # obtained via 'gene()'
    #
    # For discrete (factor) data, 'reorder' and 'relabels' are used to reorder
    # and/or rename the factor levels within the data so that the data are
    # renamed and used in the proper order by ggplot.
    #
    # *If 'mult = TRUE', takes in a list of 'target' and 'name' pairs.
    # Each 'target' must be the name of a gene or metadata.
    # Each pair is added to the dataframe with 'assay', 'slot', and
    # 'adjustment' used the same was as in the 'multi = FALSE' mode.

    if (is.null(df)) {
        df <- data.frame(row.names = .all_cells(object))
    }

    if (mult) {
        for (i in seq_along(target)) {
            df <- .add_by_cell(df, target[i], name[i], object, assay, slot,
                adjustment, reorder = NULL, relabels = NULL, mult = FALSE)
        }
    } else if (!is.null(target)) {
        # Obtain and reorder values
        values <- .var_OR_get_meta_or_gene(
            target, object, assay, slot, adjustment)
        values <- .rename_and_or_reorder(values, reorder, relabels)
        # Add
        df <- cbind(df, values)
        # Set the name
        names(df)[ncol(df)] <- name
    }

    df
}

.extract_Reduced_Dim <- function(reduction.use, dim=1, object) {
    # Extracts loadings ("embeddings") and suggested plotting label ("name")
    # for an individual dimensionality reduction dimension.

    if (is.data.frame(reduction.use) || is.matrix(reduction.use)) {
        
        embeds <- reduction.use
        
    } else if (is(object,"seurat")) {
        
        embeds <- eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@cell.embeddings")))
        colnames(embeds) <- paste0(eval(expr = parse(text = paste0(
            "object@dr$",reduction.use,"@key"))),seq_len(ncol(embeds)))
        
    } else if (is(object,"Seurat")) {
        
        .error_if_no_Seurat()
        embeds <- Seurat::Embeddings(object, reduction = reduction.use)
        
    } else if (is(object,"SingleCellExperiment")) {
        
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
    # Generates labels for plotting axes titles when use of 'reduction.use'
    # itself would not be ideal

    key <- reduction.use
    if (grepl("pca", tolower(reduction.use))) {
        key <- "PC"
    }
    if (grepl("cca", tolower(reduction.use))) {
        key <- "CC"
    }
    if (grepl("cca.aligned", tolower(reduction.use))) {
        key <- "aligned.CC"
    }
    if (grepl("ica", tolower(reduction.use))) {
        key <- "IC"
    }
    if (grepl("tsne", tolower(reduction.use))) {
        key <- "tSNE_"
    }
    key
}
