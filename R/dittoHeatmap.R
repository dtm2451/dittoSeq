########## dittoHeatmap: Builds a heatmap with given genes using pheatmap. ##########
#' Outputs a heatmap of given genes
#' @importFrom grDevices colorRampPalette
#'
#' @param genes String vector, c("gene1","gene2","gene3",...) = the list of genes to put in the heatmap. REQUIRED.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param data.type String. Options are "normalized" (data slot, default), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects). Note: scaling is performed on the data matrix by default.
#' @param heatmap.colors the colors to use within the heatmap.
#' Default is a ramp from navy to white to red with 50 slices.
#' @param scaled.to.max Logical which sets whether expression shoud be scaled between [0, 1].
#' This is recommended for single-cell datasets as they are generally enriched in 0s.
#' @param heatmap.colors.max.scaled the colors to use within the heatmap when \code{scaled.to.max} is set to \code{TRUE}.
#' Default is a ramp from white to red with 25 slices.
#' @param main String that sets the title for the heatmap.
#' @param cell.names.meta quoted "name" of a meta.data slot to use for naming the columns instead of using the raw cell/sample names.
#' @param col.annotation.metas String name of a metadata slot containing how the cells/samples should be annotated.
#' @param annotation.colors String (color) vector where each color will be assigned to an individual annotation in the generated annotation bars.
#' @param show.rownames,show.colnames Logical which sets whether rownames or colnames will be shown.
#' Note: if gene names are provided to \code{highlight.genes}, the \code{show.colnames} parameter is ignored.
#' @param data.out Logical that changes the output of the function.
#' If set to \code{TRUE}, the output will be a list containing the data that would have been used for generating the heatmap,
#' and a String showing how \code{\link[pheatmap]{pheatmap}} would have been called.
#' @param highlight.genes String vector of genes whose names you would like to show. Only these genes will be named in the resulting heatmap.
#' @param ... other arguments passed to \code{pheatmap}.
#' @description Given a set of genes, cells/samples, and metadata names for column annotations, it will retrieve the expression data for those genes and cells, and the annotation data for those cells.
#' It will then utilize these data to make a heatmap using the \code{\link[pheatmap]{pheatmap}} function of the \code{pheatmap} package.
#'
#' @note For scRNAseq data, running with \code{cluster_cols=FALSE} first is recommended because clustering of thousands of cells can tak a long time!
#'
#' @return A \code{pheatmap} object.
#' Alternatively, if \code{data.out} was set to \code{TRUE}, a list containing
#' \code{args} = a list of arguments passed to \code{pheatmap}, and
#' and \code{call} = String showing how \code{pheatmap} would have been called.
#' @seealso \code{\link[pheatmap]{pheatmap}}, for how to add additional heatmap tweaks.
#' Some examples of useful \code{pheatmap} parameters are:
#' \itemize{
#' \item \code{cluster_cols} and \code{cluster_rows} for controlling clustering
#' \item \code{treeheight_row} and \code{treeheight_col} for setting how large the trees on the side/top should be drawn.
#' \item \code{cutree_col} and \code{cutree_row} for spliting the heatmap based on kmeans clustering
#' }
#'
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = pbmc,
#'     col.annotation.metas = "ident")
#'
#' #' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     col.annotation.metas = "ident")
#'
#' # For real single cell data, you will have more cells than
#' #   in this truncated dataset,
#' #   so turning off cell clustering when trying out tweaks is recommended.
#' #   Do so by adding cluster_cols=FALSE
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = pbmc,
#'     col.annotation.metas = "ident",
#'     cluster_cols=FALSE)
#'
#' # When there are many cells, showing names becomes less useful.
#' #   Names can be turned off with the show.colnames parameter.
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = pbmc,
#'     col.annotation.metas = "ident",
#'     cluster_cols=FALSE,
#'     show.colnames = FALSE)
#'
#' # Additionally, it is recommended for single-cell data that the parameter
#' #   scaled.to.max be set to TRUE, or scaling be turned off altogether,
#' #   because these data are generally enriched for zeros that otherwise get
#' #   scaled to a negative value.
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = pbmc,
#'     col.annotation.metas = "ident",
#'     show.colnames = FALSE,
#'     scaled.to.max = TRUE)
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = pbmc,
#'     col.annotation.metas = "ident",
#'     show.colnames = FALSE,
#'     scaled.to.max = FALSE,
#'     scale = "none",
#'     heatmap.colors = colorRampPalette(c("white", "red"))(25))
#'
#' @export

dittoHeatmap <- function(
    genes=NULL, object = DEFAULT, cells.use = NULL, main = NA,
    cell.names.meta = NULL, data.type = "normalized",
    heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
    scaled.to.max = FALSE,
    heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
    col.annotation.metas = NULL,
    annotation.colors = rep(dittoColors(),10),
    data.out=FALSE, highlight.genes = NULL, show.colnames = TRUE,
    show.rownames = TRUE, ...) {

    if (is.null(genes)) {
        stop('This function is not set up to select which genes to use.\nPlease provide a set of genes.')
    }
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    #If cells.use given as logical, populate as names.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    if (is.null(main)) {
        main <- NA
    }

    #Make the data matrix
    data <- as.matrix(.which_data(data.type,object)[genes,cells.use])
    if (scaled.to.max) {
        maxs <- apply(data,1,max)
        data <- data/maxs
        heatmap.colors <- heatmap.colors.max.scaled
    }

    args <- list(
        mat = data, main = main, show_colnames = show.colnames,
        show_rownames = show.rownames, color = heatmap.colors,
        ...)

    if (is.null(args$scale)) {
        args$scale <- ifelse(scaled.to.max, "none", "row")
    }

    if (is.null(args$border_color)) {
      args$border_color <- NA
    }

    #Add cell/sample/row names unless provided separately by user
    if(is.null(args$labels_col) && !(is.null(cell.names.meta))) {
        args$labels_col <- as.character(meta(cell.names.meta, object)[all.cells %in% cells.use])
    }

    #Make a labels_row input for displaying only certain genes if genes were given to highlight.genes input
    if (!(is.null(highlight.genes)) && sum(highlight.genes %in% genes)>0) {
        highlight.genes <- highlight.genes[highlight.genes %in% genes]
        args$labels_row <- rownames(data)
        #Overwrite all non-highlight genes rownames to ""
        args$labels_row[-(match(highlight.genes,rownames(data)))] <- ""
        args$show_colnames <- TRUE
    }

    #Make the columns annotations data for colored annotation bars
    if (!is.null(col.annotation.metas)) {
        args$annotation_col <- data.frame(row.names = cells.use)
        for (i in seq_along(col.annotation.metas)) {
            args$annotation_col[,i] <- meta(col.annotation.metas[i], object)
        }
        names(args$annotation_col) <- col.annotation.metas
    }

    #Make the annotation colors
    if (is.null(args$annotation_colors)) {
        args <- .make_heatmap_annotation_colors(args, annotation.colors)
    }

    if (data.out) {
        OUT <- list(args = args, call = "do.call(pheatmap, args)")
    } else {
        OUT <- do.call(pheatmap::pheatmap, args)
    }
    OUT
}

#This next function creates pheatmap annotations_colors dataframe
    # list of character vectors, all named.
        # vector names = annotation titles
        # vector members' (colors') names = annotation identities
.make_heatmap_annotation_colors <- function(args, annotation.colors) {

    # Extract a default color-set
    annotation.colors.d <- annotation.colors
    annotation.colors.n <- rev(annotation.colors)[-seq_len(9)]

    # Initiate variables
    next.color.index.discrete <- 1
    next.color.index.numeric <- 1
    col_colors <- NULL
    row_colors <- NULL

    # Columns First (if there)
    if (!is.null(args$annotation_col)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_col,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        col_colors <- dfcolors_out$df_colors
        next.color.index.discrete <- dfcolors_out$next.color.index.discrete
        next.color.index.numeric <- dfcolors_out$next.color.index.numeric
    }

    # Rows Second (if there)
    if (!is.null(args$annotation_row)) {
        dfcolors_out <- .pick_colors_for_df(
            args$annotation_row,
            next.color.index.discrete, next.color.index.numeric,
            annotation.colors.d, annotation.colors.n)
        row_colors <- dfcolors_out$df_colors
    }

    args$annotation_colors <- c(col_colors, row_colors)
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

        # Determine the distinct contents of the first annotation
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
