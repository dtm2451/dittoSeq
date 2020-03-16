########## dittoHeatmap: Builds a heatmap with given genes using pheatmap. ##########
#' Outputs a heatmap of given genes
#' @importFrom grDevices colorRampPalette
#'
#' @param object A Seurat or SingleCellExperiment object to work with
#' @param genes String vector, c("gene1","gene2","gene3",...) = the list of genes to put in the heatmap.
#' If not provided, defaults to all genes of the object / assay.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param assay,slot single strings or integer that set which expression data to use. See \code{\link{gene}} for more information about how defaults for these are filled in when not provided.
#' @param order.by Single string or numeric vector which sets the ordering of cells/samples.
#' Can be the name of a gene, or metadata slot.
#' Alternatively, can be a numeric vector of length equal to the total number of cells/samples in object.
#' @param heatmap.colors the colors to use within the heatmap.
#' Default is a ramp from navy to white to red with 50 slices.
#' @param scaled.to.max Logical which sets whether expression shoud be scaled between [0, 1].
#' This is recommended for single-cell datasets as they are generally enriched in 0s.
#' @param heatmap.colors.max.scaled the colors to use within the heatmap when \code{scaled.to.max} is set to \code{TRUE}.
#' Default is a ramp from white to red with 25 slices.
#' @param main String that sets the title for the heatmap.
#' @param cell.names.meta quoted "name" of a meta.data slot to use for naming the columns instead of using the raw cell/sample names.
#' @param annotation.metas String name of a metadata slot containing how the cells/samples should be annotated.
#' @param annotation.colors String (color) vector where each color will be assigned to an individual annotation in the generated annotation bars.
#' @param data.out Logical that changes the output of the function.
#' If set to \code{TRUE}, the output will be a list containing the data that would have been used for generating the heatmap,
#' and a String showing how \code{pheatmap} would have been called.
#' @param highlight.genes String vector of genes whose names you would like to show. Only these genes will be named in the resulting heatmap.
#' @param cluster_cols,border_color,... other arguments passed to \code{pheatmap}.
#' @param show_colnames,show_rownames,scale,annotation_col,annotation_colors arguments passed to \code{pheatmap} that are over-ruled by certain \code{dittoHeatmap} functionality:
#' \itemize{
#' \item show_colnames (& labels_col) : if \code{cell.names.meta} is provided, pheatmaps's \code{labels_col} is utilized to show these names and \code{show_colnames} parameter is set to \code{TRUE}.
#' \item show_rownames (& labels_row) : if feature names are provided to \code{highlight.genes}, pheatmap's \code{labels_row} is utilized to show just these features' names and \code{show_rownames} parameter is set to \code{TRUE}.
#' \item scale : when parameter \code{scaled.to.max} is set to true, pheatmap's \code{scale} is set to \code{"none"} and the max scaling is performed prior to the pheatmap call.
#' \item annotation_col : when metadata names are provided \code{annotation.metas}, these are added either a barren annotation_col data.frame, or one provided by the user.
#' \item annotation_colors : dittoHeatmap fills this complicated-to-produce input in automatically by pulling from the colors given to \code{annotation.colors},
#' but it is possible to set all or some manually. dittoSeq will just fill in for any left out annotations. Format is a named (annotation_col & annotation_row colnames) list where individual values can also be named as well.
#' }
#' @description Given a set of genes, cells/samples, and metadata names for column annotations, it will retrieve the expression data for those genes and cells, and the annotation data for those cells.
#' It will then utilize these data to make a heatmap using the \code{\link[pheatmap]{pheatmap}} function of the \code{pheatmap} package.
#'
#' @return A \code{pheatmap} object.
#' Alternatively, if \code{data.out} was set to \code{TRUE}, a list containing
#' \code{args}, a list of arguments passed to \code{pheatmap}, and
#' and \code{call}, a string showing how \code{pheatmap} would have been called.
#'
#' @details
#' This function serves as a wrapper for creating pheatmap heatmaps from bulk or single-cell RNAseq data,
#' by essentially automating the data extraction and annotation building steps.
#' In order to use this function, you will need to have it installed.
#'
#' Provided with a set of \code{genes}, and an optional set of cells / samples to use (with the \code{cells.use} input),
#' the function will extract the expression matrix for each of those genes for each of those cells / samples.
#' This matrix is either left as is, default (for scaling within the ultimate call to pheatmap), or scaled by dividing each row by its maximum value (\code{scaled.to.max = TRUE}).
#'
#' When provided with a set of metadata slot names to use for building annotations (with the \code{annotation.metas} input),
#' the relevant metadata is retrieved from the \code{object} and \code{pheatmap}-ready \code{annotation_col} input is generated.
#' The input \code{annotation.colors} is used to establish the set of colors that should be used for building a \code{pheatmap}-ready \code{annotation_colors} input as well,
#' unless such an input has been provided by the user.
#' Note: Users can also provide an \code{annotation_row} dataframe for adding gene-annotations (See \code{\link[pheatmap]{pheatmap}} for details).
#' Colors for row annotations will still come from the \code{annotation.colors} input unless a \code{pheatmap}-ready \code{annotation_colors} input is provided by the user.
#'
#' When \code{data.out} is set to \code{TRUE}, a list is output instead of running pheatmap and providing the heatmap plot as output.
#' This list will contain a slot named \code{args} which includes the data matrix and all other arguments that would have been provided to the pheatmap call,
#' as well as a slot named \code{call} which shows how pheatmap would have been called.
#'
#' Many additional characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item The cells can be ordered by in a set way using the \code{order.by} input.
#' Note: It can take a long time to cluster thousands of samples or cells, so adding \code{order.by = 'useful-method'} can be quite useful.
#' \item A plot title can be added with \code{main}.
#' \item Gene or cell/sample annotations can be turned off with \code{show_rownames} and \code{show_colnames}, respectively, or...
#' \itemize{
#' \item Particular genes can also be selected for labeling using the \code{highlight.genes} input.
#' \item Names of all cells/samples replaced with the contents of a metadata slot using the \code{cell.names.meta} input.
#' }
#' \item Additional tweaks are possible through use of \code{\link[pheatmap]{pheatmap}} inputs which will be directly passed through.
#' Some examples of useful \code{pheatmap} parameters are:
#' \itemize{
#' \item \code{cluster_cols} and \code{cluster_rows} for controlling clustering.
#' Note: cluster_cols will always be over-written to be \code{FALSE} when the input \code{order.by} is used above.
#' \item \code{treeheight_row} and \code{treeheight_col} for setting how large the trees on the side/top should be drawn.
#' \item \code{cutree_col} and \code{cutree_row} for spliting the heatmap based on kmeans clustering
#' }
#' }
#'
#' @seealso
#' \code{\link[pheatmap]{pheatmap}}, for how to add additional heatmap tweaks.
#'
#' \code{\link{metaLevels}} for helping to create manual annotation_colors inputs.
#' This function universally checks the options/levels of a string, factor (filled only by default), or numerical metadata.
#'
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' # Pick a set of genes
#' genes <- getGenes(myRNA)[1:30]
#'
#' # Make a heatmap with cells/samples annotated by their clusters
#' dittoHeatmap(myRNA, genes,
#'     annotation.metas = "clustering")
#'
#' # For single-cell data, you will typically have more cells than can be
#' # clustered quickly. Thus, cell clustering is turned off by default for
#' # single-cell data.
#'
#' # Using the 'order.by' input:
#' #   ordering by a useful metadata or gene is generally more helpful
#' #   For single-cell data, order.by defaults to the first element given to
#' #     annotation.metas.
#' #   For bulk data, order.by must be set separately.
#' dittoHeatmap(myRNA, genes,
#'     annotation.metas = "clustering",
#'     order.by = "clustering")
#'
#' # When there are many cells, showing names becomes less useful.
#' #   Names can be turned off with the show_colnames parameter.
#' dittoHeatmap(myRNA, genes,
#'     annotation.metas = "groups",
#'     order.by = "groups",
#'     show_colnames = FALSE)
#'
#' # Additionally, it is recommended for single-cell data that the parameter
#' #   scaled.to.max be set to TRUE, or scale be "none" and turned off altogether,
#' #   because these data are generally enriched for zeros that otherwise get
#' #   scaled to a negative value.
#' dittoHeatmap(myRNA, genes, annotation.metas = "groups",
#'     order.by = "groups", show_colnames = FALSE,
#'     scaled.to.max = TRUE)
#' dittoHeatmap(myRNA, genes, annotation.metas = "groups",
#'     order.by = "groups", show_colnames = FALSE,
#'     scaled.to.max = FALSE,
#'     heatmap.colors = colorRampPalette(c("white", "red"))(25))
#'
#' @author Daniel Bunis
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

dittoHeatmap <- function(
    object, genes=getGenes(object), cells.use = NULL,
    cluster_cols = .is_bulk(object),
    annotation.metas = NULL,
    order.by = .default_order(object, annotation.metas),
    main = NA, cell.names.meta = NULL,
    assay = .default_assay(object), slot = .default_slot(object),
    heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
    scaled.to.max = FALSE,
    heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
    annotation.colors = c(rep(dittoColors(),9),dittoColors()[seq_len(7)]),
    annotation_col = NULL, annotation_colors = NULL,
    data.out=FALSE, highlight.genes = NULL, show_colnames = TRUE,
    show_rownames = TRUE, scale = "row", border_color = NA, ...) {

    # If cells.use given as logical, populate as names.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    # Adjust title (off = NA in pheatmap)
    if (is.null(main)) {
        main <- NA
    }

    # Make the data matrix
    data <- as.matrix(.which_data(assay,slot,object)[genes,cells.use])
    if (sum(rowSums(data)==0)>0) {
        data <- data[rowSums(data)!=0,]
        if (nrow(data)==0) {
            stop("No target genes are expressed in the 'cells.use' subset")
        }
        warning("Gene(s) removed due to absence of expression within the 'cells.use' subset")
    }

    # Make the columns annotations data for colored annotation bars
    if (is.null(annotation_col)) {
        annotation_col <- data.frame(row.names = cells.use)
    } else if (!all(cells.use %in% rownames(annotation_col))) {
        stop("rows of 'annotation_col' must be cell/sample names of all cells/samples being displayed")
    }
    if (!is.null(annotation.metas)) {
        annotation_col <- rbind(
            as.data.frame(getMetas(object, names.only = FALSE)[cells.use, annotation.metas, drop = FALSE]),
            annotation_col[cells.use, , drop=FALSE])
    }

    # Create the base pheatmap inputs
    args <- list(
        mat = data, main = main, show_colnames = show_colnames,
        show_rownames = show_rownames, color = heatmap.colors,
        cluster_cols = cluster_cols, border_color = border_color,
        scale = scale, ...)
    # Add annotation_col / annotation_colors only if needed
    if (ncol(annotation_col)>0 || !is.null(args$annotation_row)) {
        if (ncol(annotation_col)>0) {
            args$annotation_col <- annotation_col
        }
        args$annotation_colors <- annotation_colors
        args <- .make_heatmap_annotation_colors(args, annotation.colors)
    }

    if (scaled.to.max) {
        maxs <- apply(args$mat,1,max)
        args$mat <- args$mat/maxs
        args$color <- heatmap.colors.max.scaled
        args$scale <- "none"
    }
    if (!is.null(order.by)) {
        order_data <- .var_OR_get_meta_or_gene(order.by, object, assay, slot)
        args$mat <- args$mat[,order(order_data[all.cells %in% cells.use])]
    }
    # Make a labels_row input for displaying only certain genes if genes were given to 'highlight.genes'
    if (!(is.null(highlight.genes)) && sum(highlight.genes %in% rownames(data))>0) {
        highlight.genes <- highlight.genes[highlight.genes %in% rownames(data)]
        args$labels_row <- rownames(data)
        #Overwrite all non-highlight genes rownames to ""
        args$labels_row[-(match(highlight.genes,rownames(data)))] <- ""
        args$show_rownames <- TRUE
    }

    # Add cell/sample/row names unless provided separately by user
    if(is.null(args$labels_col) && !(is.null(cell.names.meta))) {
        names <- .var_OR_get_meta_or_gene(cell.names.meta, object)
        args$labels_col <- as.character(names[colnames(args$mat)])
    }

    if (data.out) {
        OUT <- list(args = args, call = "do.call(pheatmap, args)")
    } else {
        OUT <- do.call(pheatmap::pheatmap, args)
    }
    OUT
}

# This next function creates pheatmap annotations_colors dataframe
    # list of character vectors, all named.
        # vector names = annotation titles
        # vector members' (colors') names = annotation identities
.make_heatmap_annotation_colors <- function(args, annotation.colors) {

    # Extract a default color-set
    annotation.colors.d <- annotation.colors
    annotation.colors.n <- rev(annotation.colors)

    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    col_colors <- NULL
    row_colors <- NULL
    user.provided <- names(args$annotation_colors)

    # Columns First (if there)
    if (!is.null(args$annotation_col) && ncol(args$annotation_col)>0) {
        make.these <- !(colnames(args$annotation_col) %in% user.provided)
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_col[,make.these, drop = FALSE],
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        col_colors <- dfcolors_out$df_colors
        next.color.index.discrete <- dfcolors_out$next.color.index.discrete
        next.color.index.numeric <- dfcolors_out$next.color.index.numeric
    }

    # Rows Second (if there)
    if (!is.null(args$annotation_row)) {
        make.these <- !(colnames(args$annotation_row) %in% user.provided)
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_row[,make.these, drop = FALSE],
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        row_colors <- dfcolors_out$df_colors
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
    annotation.colors.d, annotation.colors.n
    ) {

    df_colors <- NULL
    for (i in seq_len(ncol(annotation_df))){

        # Determine the distinct contents of the annotation
        in.this.annot <- levels(as.factor(annotation_df[,i]))

        # Make new colors
        if(!is.numeric(annotation_df[,i])){
            # Take colors for each, and name them.
            new.colors <- annotation.colors.d[
                seq_along(in.this.annot) + next.color.index.discrete - 1
                ]
            names(new.colors) <- in.this.annot

            next.color.index.discrete <-
                next.color.index.discrete + length(in.this.annot)
        } else {
            # Make a 100 color split as in pheatmap code.
            a <- cut(
                annotation_df[order(annotation_df[,i]),i],
                breaks = 100)
            # Assign to colors.
            this.ramp <- annotation.colors.n[next.color.index.numeric]
            new.colors <-
                grDevices::colorRampPalette(c("white",this.ramp))(100)[a]

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

.default_order <- function(object, annotation.metas) {
    if (!is.null(annotation.metas) && !.is_bulk(object)) {
        return(annotation.metas[1])
    } else {
        return(NULL)
    }
}

