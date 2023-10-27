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
        stop("No genes or metadata names added to 'hover.data'")
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
    # Creates a single character vector where each element is the hoverstring
    # for a given row of the provided 'df' with structure
    # "col1name: col1-value\ncol2name: var2-value\ncol3name: var3-value\n..."
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
    # Takes in string vector or factor 'orig.data', integer vector 'reorder',
    # and string vector 'relabels'.
    # Turns character vectors into factors
    # Reorders the level of the factor based on indices provided to 'reorder'
    # Re-labels the levels of the factor based on lebels provided to 'relabels'
    
    # Need to skip if no reorder or relabels demands because without them, the factor() call at the end will drop unused levels!
    if (is.numeric(orig.data) || is.null(reorder) && is.null(relabels)) {
        return(orig.data)
    }
    rename.args <- list(x = orig.data)
    if (!(is.null(reorder))) {
        if (length(reorder)!=length(levels(factor(orig.data)))) {
            stop("incorrect number of indices provided to 'reorder' input")
        }
        rename.args$levels <- levels(factor(orig.data))[reorder]
    }
    if (!(is.null(relabels))) {
        if (length(relabels)!=length(levels(factor(orig.data)))) {
            stop("incorrect number of labels provided to 'relabel' input")
        }
        rename.args$labels <- relabels
    }
    do.call(factor, args = rename.args)
}

.keep_levels_if_factor <- function(current.values, orig.structure, drop.unused = TRUE) {
    # current.values should be the column as it exists in the to-plot-dataframe
    # orig.structure should be the original values from the object, already subset by cells.use
    if (is.factor(orig.structure)) {
        levs_use <- levels(orig.structure)
        if (drop.unused) {
            levs_keep <- unique(as.character(orig.structure))
            levs_use <- levs_use[levs_use %in% levs_keep]
        }
        current.values <- factor(
            current.values,
            levels = levs_use
        )
    }
    current.values
}

.swap_rownames <- function(object, swap.rownames = NULL) {

    if (identical(swap.rownames, NULL)) {
        return(object)
    }

    # SummarizedExperiment / SCEs only
    if (is(object, "SummarizedExperiment") ) {

        # Alternate experiments, simplified usage = find unnamed swap.rownames elements, prioritized by order, in all experiments and swap to them.
        if (is(object,"SingleCellExperiment") && length(SingleCellExperiment::altExpNames(object))>0 && identical(NULL, names(swap.rownames))) {

            .which_swap_rownames <- function(swaps = swap.rownames, exp, name) {
                use <- swaps[swaps %in% names(SummarizedExperiment::rowData(exp))]
                if (length(use)==0) {
                    warning("Skipping swapping rownames for ", name," experiment because no 'swap.rownames' elements are a rowData column name.")
                    return(NULL)
                } else {
                    use[1]
                }
            }
            # Determine explicit swaps needed for main and alternate experiments
            # main
            new_swaps <- c('main' = .which_swap_rownames(swap.rownames, object, 'main'))
            # alternates
            for (expName in SingleCellExperiment::altExpNames(object)) {
                next_new <- .which_swap_rownames(swap.rownames, SingleCellExperiment::altExp(object, expName), expName)
                if (!identical(next_new, NULL)) {
                    names(next_new) <- expName
                    new_swaps <- c(new_swaps, next_new)
                }
            }

            # Execute with determined explicit swaps
            return(.swap_rownames(object, new_swaps))
        }

        # Recursion if length > 1
        if (length(swap.rownames) > 1) {
            for (i in seq_along(swap.rownames)) {
                object <- .swap_rownames(object, swap.rownames[i])
            }
            return(object)
        }

        # alternate experiments, explicit usage
        if (!identical(NULL, names(swap.rownames)) && !names(swap.rownames) %in% c('', 'main') && is(object,"SingleCellExperiment")) {
            if (!is(object,"SingleCellExperiment")) {
                warning("Ignoring names in 'swap.rownames' for non-SCE 'object'.")
            } else {
                if (names(swap.rownames) == 'altexp') {
                    names(swap.rownames) <- SingleCellExperiment::altExpNames(object)[1]
                }
                if (!names(swap.rownames) %in% SingleCellExperiment::altExpNames(object)) {
                    stop("A 'swap.rownames' element has a name that is not either the name of an alternate experiment within 'object' or the generic \"altexp\".")
                }
                if (!swap.rownames %in% colnames(SummarizedExperiment::rowData(SingleCellExperiment::altExp(object, names(swap.rownames))))) {
                    stop("A 'swap.rownames' value is not a rowData column of this object's ", names(swap.rownames), " alternate experiment.")
                }
                rownames(SingleCellExperiment::altExp(object, names(swap.rownames))) <- rowData(SingleCellExperiment::altExp(object, names(swap.rownames)))[,swap.rownames]
                return(object)
            }
        }

        # 'main' experiments and SEs
        if (!swap.rownames %in% names(SummarizedExperiment::rowData(object))) {
            stop("'swap.rownames' is not a column of 'rowData(object)'")
        }
        rownames(object) <- rowData(object)[,swap.rownames]

    }

    object
}


