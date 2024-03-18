########## dittoHeatmap: Builds a heatmap with given genes using pheatmap. ##########
#' Outputs a heatmap of given genes
#' @importFrom grDevices colorRampPalette
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param genes String vector, c("gene1","gene2","gene3",...) = the list of genes to put in the heatmap.
#' If not provided, defaults to all genes of the object / assay.
#' @param metas String vector, c("meta1","meta2","meta3",...) = the list of metadata variables to put in the heatmap.
#' @param complex Logical which sets whether the heatmap should be generated with ComplexHeatmap (\code{TRUE}) versus pheatmap (\code{FALSE}, default).
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param order.by Single string, string vector, or numeric vector which sets how cells/samples (columns) will be ordered when \code{cluster_cols = FALSE}.
#' 
#' Strings should be the name of a gene, or metadata slot, but can also be multiple such values in order of priority.
#' 
#' Alternatively, can be a numeric vector which gives the column index order directly.
#' @param heatmap.colors the colors to use within the heatmap when (default setting) \code{scaled.to.max} is set to \code{FALSE}.
#' Default is a ramp from navy to white to red with 50 slices.
#' @param scaled.to.max Logical, \code{FALSE} by default, which sets whether expression shoud be scaled between [0, 1].
#' This is recommended for single-cell datasets as they are generally enriched in 0s.
#' @param heatmap.colors.max.scaled the colors to use within the heatmap when \code{scaled.to.max} is set to \code{TRUE}.
#' Default is a ramp from white to red with 25 slices.
#' @param main String that sets the title for the heatmap.
#' @param cell.names.meta quoted "name" of a meta.data slot to use for naming the columns instead of using the raw cell/sample names.
#' @param annot.by String name of any metadata slots containing how the cells/samples should be annotated.
#' @param annot.colors String (color) vector where each color will be assigned to an individual annotation in the generated annotation bars.
#' @param data.out Logical. When set to \code{TRUE}, changes the output from the heatmat itself, to a list containing all arguments that would have be passed to \code{\link{pheatmap}} for heatmap generation.
#' (Can be useful for troubleshooting or customization.)
#' @param highlight.features String vector of genes/metadata whose names you would like to show. Only these genes/metadata will be named in the resulting heatmap.
#' @param cluster_cols,border_color,legend_breaks,breaks,drop_levels,... other arguments passed to \code{\link[pheatmap]{pheatmap}} directly (or to \code{\link[ComplexHeatmap]{pheatmap}} if \code{complex = TRUE}).
#' @param show_colnames,show_rownames,scale,annotation_col,annotation_colors arguments passed to \code{pheatmap} that are over-ruled by certain \code{dittoHeatmap} functionality:
#' \itemize{
#' \item show_colnames (& labels_col): if \code{cell.names.meta} is provided, pheatmaps's \code{labels_col} is utilized to show these names and \code{show_colnames} parameter is set to \code{TRUE}.
#' \item show_rownames (& labels_row): if feature names are provided to \code{highlight.features}, pheatmap's \code{labels_row} is utilized to show just these features' names and \code{show_rownames} parameter is set to \code{TRUE}.
#' \item scale: when parameter \code{scaled.to.max} is set to true, pheatmap's \code{scale} is set to \code{"none"} and the max scaling is performed prior to the pheatmap call.
#' \item annotation_col: Can be provided as normal by the user and any metadata given to \code{annot.by} will then be appended.
#' \item annotation_colors: dittoHeatmap fills this complicated-to-produce input in automatically by pulling from the colors given to \code{annot.colors},
#' but it is possible to set all or some manually. dittoSeq will just fill any left out annotations. Format is a named (annotation_col & annotation_row colnames) character vector list where individual color values can also be named.
#' }
#'
#' @inheritParams dittoPlot
#'
#' @description Given a set of genes, cells/samples, and metadata names for column annotations, this function will retrieve the expression data for those genes and cells, and the annotation data for those cells.
#' It will then utilize these data to make a heatmap using the \code{\link[pheatmap]{pheatmap}} function of either the \code{pheatmap} (default) or \code{ComplexHeatmap} package.
#'
#' @return A \code{pheatmap} object.
#'
#' Alternatively, if \code{complex} is set to \code{TRUE}, a \code{\link[ComplexHeatmap]{Heatmap}}
#' 
#' Alternatively, if \code{data.out} is set to \code{TRUE}, a list containing all arguments that would have be passed to pheatmap to generate such a heatmap.
#'
#' @details
#' This function serves as a wrapper for creating heatmaps from bulk or single-cell RNAseq data with pheatmap::\code{\link{pheatmap}},
#' by essentially automating the data extraction and annotation building steps.
#' (Or alternatively with ComplexHeatmap::\code{\link[ComplexHeatmap]{pheatmap}} if \code{complex} is set to \code{true}.
#'
#' The function will extract the expression matrix for a set of \code{genes} and/or an optional subset of cells / samples to use via \code{cells.use},
#' This matrix is either left as is, default (for scaling within the ultimate call to pheatmap), or if \code{scaled.to.max = TRUE}, is scaled by dividing each row by its maximum value.
#'
#' When provided with a set of metadata slot names to use for building annotations (with the \code{annot.by} input),
#' the relevant metadata is retrieved from the \code{object} and compiled into a \code{pheatmap}-ready \code{annotation_col} input.
#' The input \code{annot.colors} is used to establish the set of colors that should be used for building a \code{pheatmap}-ready \code{annotation_colors} input as well,
#' unless such an input has been provided by the user. See below for further details.
#'
#' @section Many additional characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item The cells can be ordered in a set way using the \strong{\code{order.by}} input.
#'
#' Such ordering happens by default for single-cell RNAseq data when any metadata are provided to \code{annot.by} as it is often unfeasible to cluster thousands of cells.
#' \item A plot title can be added with \strong{\code{main}}.
#' \item Gene or cell/sample names can be hidden with \code{show_rownames} and \code{show_colnames}, respectively, or...
#' \itemize{
#' \item Particular features can also be selected for labeling using the \code{highlight.features} input.
#' \item Names of all cells/samples can be replaced with the contents of a metadata slot using the \code{cell.names.meta} input.
#' }
#' \item Additional tweaks are possible through use of \code{\link[pheatmap]{pheatmap}} inputs which will be directly passed through.
#' Some examples of useful \code{pheatmap} parameters are:
#' \itemize{
#' \item \code{cluster_cols} and \code{cluster_rows} for controlling clustering.
#' Note: cluster_cols will always be over-written to be \code{FALSE} when the input \code{order.by} is used above.
#' \item \code{treeheight_row} and \code{treeheight_col} for setting how large the trees on the side/top should be drawn.
#' \item \code{cutree_col} and \code{cutree_row} for spliting the heatmap based on kmeans clustering
#' }
#' \item When \code{complex} is set to \code{TRUE}, additional inputs for the \code{\link[ComplexHeatmap]{Heatmap}} function can be given as well.
#' Some examples:
#' \itemize{
#' \item \code{use_raster} to have the heatmap rasterized/flattened to pixels which can make working with large heatmaps in a figure editor, like Illustrator, simpler.
#' \item \code{name} to give the heatmap color scale a custom title.
#' }
#' }
#'
#' @section Customized annotations:
#' In typical operation, dittoHeatmap pulls metadata annotations given to \code{annot.by} to build a pheatmap-\code{annotation_col} input,
#' then it uses the colors provided to \code{annot.colors} to create the pheatmap-\code{annotation_colors} input which sets the annotation coloring.
#' Specifically...
#' \itemize{
#' \item colors for the values of \strong{discrete} metadata are pulled from the \emph{start} of the \code{annot.colors} vector, in the order that they are given to \code{annot.by}
#' \item colors for the values of \strong{continuous} metadata are pulled from the \emph{end} of the \code{annot.colors} vector, in the order that they are given to \code{annot.by}
#' }
#'
#' To customize colors or add additional column or row annotations, users can also provide \code{annotation_colors}, \code{annotation_col}, or \code{annotation_row} pheatmap-inputs directly.
#' General structure is described below, but see \code{\link[pheatmap]{pheatmap}} for additional details and examples.
#' \itemize{
#' \item \code{annotation_col} = a data.frame with rownames of the barcodes/names of all cells/samples in the dataset & columns representing annotations.
#' Names of columns are used as the annotation titles. *dittoSeq will append any \code{annot.by} annotations to this dataframe.
#' \item \code{annotation_row} = a data.frame with rownames of the genes/feature of the dataset & columns representing annotations.
#' Names of columns are used as the annotation titles.
#' \item \code{annotation_colors} = a named list of string (color) vectors.
#' Vectors must be named by the row or column annotation title that they are associated with.
#' Optionally, individual colors can be named with the values that they should be associated with.
#'
#' Partial \code{annotation_colors} lists (containing vectors for only certain annotations) will have colors for left out annotations filled in automatically.
#' For such filling, \code{annot.colors} are pulled for column annotations first, then for row annotations.
#' }
#'
#' @seealso
#' pheatmap::\code{\link[pheatmap]{pheatmap}}, for how to add additional heatmap tweaks,
#' OR or ComplexHeatmap::\code{\link[ComplexHeatmap]{pheatmap}} and \code{\link[ComplexHeatmap]{Heatmap}} for when you want to turn on rasterization or any additional customizations offered by this fantastic package.
#'
#' \code{\link{metaLevels}} for helping to create manual annotation_colors inputs.
#' This function universally checks the options/levels of a string, factor (filled only by default), or numerical metadata.
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' scRNA <- setBulk(myRNA, FALSE)
#' 
#' # We now have two SCEs for our example purposes:
#'   # 'myRNA' will be treated as a bulk RNAseq dataset
#'   # 'scRNA' will be treated as a single-cell RNAseq dataset
#'
#' # Pick a set of genes
#' genes <- getGenes(myRNA)[1:30]
#'
#' # Make a heatmap with cells/samples annotated by their clusters
#' dittoHeatmap(myRNA, genes,
#'     annot.by = "clustering")
#'
#' # For single-cell data, you will typically have more cells than can be
#' # clustered quickly. Thus, cell clustering is turned off by default for
#' # single-cell data.
#' dittoHeatmap(scRNA, genes,
#'     annot.by = "clustering")
#'
#' # Using the 'order.by' input:
#' #   Ordering by a useful metadata or gene is often helpful.
#' #   For single-cell data, order.by defaults to the first element given to
#' #     annot.by.
#' #   For bulk data, order.by must be set separately.
#' dittoHeatmap(myRNA, genes,
#'     annot.by = "clustering",
#'     order.by = "clustering",
#'     cluster_cols = FALSE)
#' # 'order.by' can be multiple metadata/genes, or a vector of indexes directly 
#' dittoHeatmap(scRNA, genes,
#'     annot.by = "clustering",
#'     order.by = c("clustering", "timepoint"))
#' dittoHeatmap(scRNA, genes,
#'     annot.by = "clustering",
#'     order.by = ncol(scRNA):1)
#'
#' # When there are many cells, showing names becomes less useful.
#' #   Names can be turned off with the 'show_colnames' parameter.
#' dittoHeatmap(scRNA, genes,
#'     annot.by = "groups",
#'     show_colnames = FALSE)
#' 
#' # When theree are many many cells & genes, rasterization can be super useful
#' # as well.
#' #   Rasterization, or flattening of the distinct color objects to a matrix of
#' #   pixels, is the default for large heatmaps in the ComplexHeatmap package,
#' #   and you can have the heatmap rendered with this package (rather than the
#' #   pheatmap package) by setting 'complex = TRUE'.
#' #   Our data here is too small to hit that defaulting switch, so lets give
#' #   the direct input, 'use_raster' as well:
#' if (requireNamespace("ComplexHeatmap")) { # Checks if you have the package.
#'     dittoHeatmap(scRNA, genes, annot.by = "groups", show_colnames = FALSE,
#'         complex = TRUE,
#'         use_raster = TRUE)
#' }
#'
#' # Additionally, it is recommended for single-cell data that the parameter
#' #   scaled.to.max be set to TRUE, or scale be "none" and turned off altogether,
#' #   because these data are generally enriched for zeros that otherwise get
#' #   scaled to a negative value.
#' dittoHeatmap(myRNA, genes, annot.by = "groups",
#'     order.by = "groups", show_colnames = FALSE,
#'     scaled.to.max = TRUE)
#' 
#'
#' @author Daniel Bunis and Jared Andrews
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

dittoHeatmap <- function(
    object,
    genes = getGenes(object, assay),
    metas = NULL,
    cells.use = NULL,
    annot.by = NULL,
    order.by = .default_order(object, annot.by),
    main = NA,
    cell.names.meta = NULL,
    assay = .default_assay(object),
    slot = .default_slot(object),
    swap.rownames = NULL,
    heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
    scaled.to.max = FALSE,
    heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
    annot.colors = c(dittoColors(),dittoColors(1)[seq_len(7)]),
    annotation_col = NULL,
    annotation_colors = NULL,
    data.out=FALSE,
    highlight.features = NULL,
    show_colnames = isBulk(object),
    show_rownames = TRUE,
    scale = "row",
    cluster_cols = isBulk(object),
    border_color = NA,
    legend_breaks = NA,
    drop_levels = FALSE,
    breaks = NA,
    complex = FALSE,
    ...) {

    # If cells.use given as logical, populate as names.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    ### Obtain all needed data
    object <- .swap_rownames(object, swap.rownames)
    
    data <- .get_heatmap_data(object, genes, metas, assay, slot, cells.use)
    
    if (!is.null(cell.names.meta)) {
        cell.names <- .var_OR_get_meta_or_gene(cell.names.meta, object)
    } else {
        cell.names <- NULL
    }
    
    if (!is.null(order.by)) {
        
        if (is.numeric(order.by)) {
            ordering <- order.by
            # Trim by cells.use if longer than cells.use
            if (length(ordering) > length(cells.use)){
                ordering <- ordering[all.cells %in% cells.use]
            }
        } else {
            order_data <- lapply(
                order.by,
                function(x) {
                    .var_OR_get_meta_or_gene(x, object, assay, slot)[all.cells %in% cells.use]
                })
            ordering <- do.call(order, order_data)
        }
        
        data <- data[, ordering]
    }
    
    # Make the columns annotations data
    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = cells.use)
    } else if (!all(cells.use %in% rownames(annotation_col))) {
        stop("rows of 'annotation_col' must be cell/sample names of all cells/samples being displayed")
    }
    if (!is.null(annot.by)) {
        if ("ident" %in% annot.by && isMeta("ident", object) && !"ident" %in% getMetas(object)) {
            object$ident <- meta("ident", object)
        }
        annotation_col <- rbind(
            as.data.frame(getMetas(object, names.only = FALSE)[cells.use, annot.by, drop = FALSE]),
            annotation_col[cells.use, , drop=FALSE])
    }
    
    # Set title (off = NA in pheatmap)
    if (is.null(main)) {
        main <- NA
    }
    
    ### Prep inputs for heatmap contructor calls
    args <- .prep_ditto_heatmap(
        data, cells.use, all.cells, cell.names, main,
        heatmap.colors, scaled.to.max, heatmap.colors.max.scaled, annot.colors,
        annotation_col, annotation_colors, highlight.features, show_colnames,
        show_rownames, scale, cluster_cols, border_color, legend_breaks,
        drop_levels, breaks, ...)
    
    if (data.out) {
        OUT <- args
    } else if (complex) {
        .error_if_no_complexHm()
        OUT <- do.call("pheatmap", args, envir = asNamespace("ComplexHeatmap"))
    } else {
        OUT <- do.call(pheatmap::pheatmap, args)
    }
    OUT
}

.prep_ditto_heatmap <- function(
    data,
    cells.use,
    all.cells,
    cell.names,
    main,
    heatmap.colors,
    scaled.to.max,
    heatmap.colors.max.scaled,
    annot.colors,
    annotation_col,
    annotation_colors,
    highlight.features,
    show_colnames,
    show_rownames,
    scale,
    cluster_cols,
    border_color,
    legend_breaks,
    drop_levels,
    breaks,
    ...) {
    
    # Create the base pheatmap inputs
    args <- list(
        mat = data, main = main, show_colnames = show_colnames,
        show_rownames = show_rownames, color = heatmap.colors,
        cluster_cols = cluster_cols, border_color = border_color,
        scale = scale, breaks = breaks, legend_breaks = legend_breaks,
        drop_levels = drop_levels, ...)
    
    # Adjust data
    
    if (scaled.to.max) {
        args <- .scale_to_max(args, heatmap.colors.max.scaled)
    }
    
    # Add annotation_col / annotation_colors only if needed
    if (ncol(annotation_col)>0 || "annotation_row" %in% names(args)) {
        if (ncol(annotation_col)>0) {
            args$annotation_col <- annotation_col[colnames(args$mat),, drop = FALSE]
        }
        args$annotation_colors <- annotation_colors
        # Add any missing annotation colors
        args <- .make_heatmap_annotation_colors(args, annot.colors, drop_levels)
    }

    # Make a labels_row input for displaying only certain genes/variables if given to 'highlight.features'
    if (!(is.null(highlight.features)) && sum(highlight.features %in% rownames(data))>0) {
        highlight.features <- highlight.features[highlight.features %in% rownames(data)]
        args$labels_row <- rownames(data)
        #Overwrite all non-highlight genes rownames to ""
        args$labels_row[-(match(highlight.features,rownames(data)))] <- ""
        args$show_rownames <- TRUE
    }

    # Add cell/sample/row names unless provided separately by user
    if(!"labels_col" %in% names(args) && !is.null(cell.names)) {
        args$labels_col <- as.character(cell.names[colnames(args$mat)])
        args$show_colnames <- TRUE
    }
    
    args
}

.scale_to_max <- function(args, heatmap.colors.max.scaled) {
    # Max scales the data matrix,
    # addss max-sclaing-colors
    # adjusts the color and legend breaks accordingly (unless set by user)
    maxs <- apply(args$mat,1,max)
    args$mat <- args$mat/maxs
    args$scale <- "none"
    args$color <- heatmap.colors.max.scaled
    if (identical(args$legend_breaks, NA)) {
        args$legend_breaks <- seq(0, 1, 0.2)
    }
    if (identical(args$breaks, NA)) {
        args$breaks <- seq(0, 1, length.out = length(args$color)+1)
    }

    args
}

# This next function creates pheatmap annotations_colors dataframe
    # list of character vectors, all named.
        # vector names = annotation titles
        # vector members' (colors') names = annotation identities
.make_heatmap_annotation_colors <- function(args, annot.colors, drop_levels) {

    # Extract a default color-set
    annot.colors.d <- annot.colors
    annot.colors.n <- rev(annot.colors)

    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    col_colors <- NULL
    row_colors <- NULL
    user.provided <- names(args$annotation_colors)

    # Columns First (if there)
    if ("annotation_col" %in% names(args) && ncol(args$annotation_col)>0) {
        make.these <- !(colnames(args$annotation_col) %in% user.provided)
        if (any(make.these)) {
            dfcolors_out <- .pick_colors_for_df(
                args$annotation_col[,make.these, drop = FALSE],
                next.color.index.discrete, next.color.index.numeric,
                annot.colors.d, annot.colors.n, drop_levels)
            col_colors <- dfcolors_out$df_colors
            next.color.index.discrete <- dfcolors_out$next.color.index.discrete
            next.color.index.numeric <- dfcolors_out$next.color.index.numeric
        }
    }

    # Rows Second (if there)
    if ("annotation_row" %in% names(args) && ncol(args$annotation_row)>0) {
        make.these <- !(colnames(args$annotation_row) %in% user.provided)
        if (any(make.these)) {
            dfcolors_out <- .pick_colors_for_df(
                args$annotation_row[,make.these, drop = FALSE],
                next.color.index.discrete, next.color.index.numeric,
                annot.colors.d, annot.colors.n, drop_levels)
            row_colors <- dfcolors_out$df_colors
        }
    }

    # Combine new with user.provided
    args$annotation_colors <- c(
        args$annotation_colors, col_colors, row_colors)

    args
}

# Interpret annotations dataframe,
# Pick, name, and add colors.
.pick_colors_for_df <- function(
    annotation_df,
    next.color.index.discrete, next.color.index.numeric,
    annot.colors.d, annot.colors.n, drop_levels
    ) {

    df_colors <- list()
    for (i in seq_len(ncol(annotation_df))){

        # Make new colors
        if(!is.numeric(annotation_df[,i])){
            # Discrete, determine the distinct contents of the annotation first
            in.this.annot <- levels(as.factor(annotation_df[,i]))
            if(drop_levels){
                levels.with.data <- unique(as.character(annotation_df[,i]))
                in.this.annot <- in.this.annot[in.this.annot %in% levels.with.data]
            }
            
            # Take colors for each, and name them.
            new.colors <- annot.colors.d[
                seq_along(in.this.annot) + next.color.index.discrete - 1
                ]
            names(new.colors) <- in.this.annot

            next.color.index.discrete <-
                next.color.index.discrete + length(in.this.annot)
        } else {
            # Numeric, just need colors for min (white) and max
            new.colors <-c("white",annot.colors.n[next.color.index.numeric])
            
            next.color.index.numeric <- next.color.index.numeric + 1
        }
        
        # Add the new colors as the list
        df_colors <- c(
            df_colors,
            list(new.colors))
    }
    names(df_colors) <- names(annotation_df)
    
    list(df_colors = df_colors,
        next.color.index.discrete = next.color.index.discrete,
        next.color.index.numeric = next.color.index.numeric)
}

# Get the heatmap data matrix.
.get_heatmap_data <- function(object, genes, metas, assay, slot, cells.use) {
    if (!is.null(genes)) { 
        data <- as.matrix(.which_data(assay,slot,object)[genes,cells.use]) 
    } else { 
        data <- NULL 
    } 

    if (!is.null(metas)) { 
        met.data <- as.matrix(t(getMetas(object, names.only = FALSE)[cells.use, metas])) 
        data <- rbind(data, met.data) 
    } 

    if (is.null(data)) { 
        stop("No 'genes' or 'metas' requested") 
    } 

    if (any(rowSums(data!=0)==0)) {
        data <- data[rowSums(data!=0)!=0,]
        if (nrow(data)==0) { 
            stop("No target genes/metadata features have non-zero values in the 'cells.use' subset") 
        } 
        warning("Gene(s) or metadata removed due to absence of non-zero values within the 'cells.use' subset") 
    } 

    data
}
