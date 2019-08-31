########## dittoHeatmap: Builds a heatmap with given genes using pheatmap. ##########
#' Outputs a heatmap of given genes
#' @importFrom grDevices colorRampPalette
#'
#' @param genes c("gene1","gene2","gene3",...) = list of genes to put in the heatmap. REQUIRED.
#' @param object String name of the object to draw from.  Alternatively, the \code{\link[Seurat]{Seurat}},  or \code{\linkS4class{RNAseq}}, or \code{\link[SingleCellExperiment]{SingleCellExperiment}} object itself. REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param cells.use String vector of cell/sample names to include, OR logical vector that is the same length as the number of cells in the object.
#' @param data.type String. Options are "normalized" (data slot, default), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects). Note: scaling is performed on the data matrix by default.
#' @param heatmap.colors the colors to use within the heatmap.
#' Default is a ramp from navy to white to red with 50 slices.
#' @param scale.to.max Logical which sets whether expression shoud be scaled between [0, 1].
#' This is recommended for single-cell datasets as they are generally enriched in 0s.
#' @param heatmap.colors.max.scaled the colors to use within the heatmap when \code{scaled.to.max} is set to \code{TRUE}.
#' Default is a ramp from white to red with 25 slices.
#' @param main String that sets the title for the heatmap.
#' @param cell.names.meta quoted "name" of a meta.data slot to use for naming the columns instead of the raw cell/sample names.
#' @param col.annotation.metas String name of a metadata slot containing how the cells/samples should be annotated. Default = \code{NULL}.
#' @param annotation.colors List of String (color) vectors where each vector is the set of colors for an individual annotation bar.
#' @param show.rownames,show.colnames Logical which sets whether rownames or colnames will be shown.
#' Note: if gene names are provided to \code{highlight.genes}, the \code{show.colnames} parameter is ignored.
#' @param data.out = If set to TRUE, the output will be a list containing \code{args},
#' a list of arguments passed to \code{pheatmap},
#' and \code{call} containing how \code{pheatmap} would have been called.
#' @param highlight.genes = A list of genes whose names you would like to show, c("Gene1","Gene2","Gene3").  Only these genes will be named in the resulting heatmap.
#' @param ... other arguments that will be passed to \code{\link[pheatmap]{pheatmap}}.
#' For scRNAseq data, I recommend running with \code{cluster_cols=FALSE} first because clustering of thousands of cells can tak a long time!
#' Common options are \code{main} for adding a title,\code{cluster_cols} or \code{cluster_rows} for controlling clustering, \code{treeheight_row} or \code{treeheight_col} for setting how large the trees on the side/top should be drawn.
#' @description This function is a wrapper for the \code{\link[pheatmap]{pheatmap}} function of the \code{pheatmap} package.
#' Given a set of genes, cells, and annotation metadata names, it will retrieve the expression data for those genes and cells, and the annotation data for those cells,
#' then it will utilize the \code{pheatmap} function to make a heatmap.
#' @details
#'
#' Notes on adding row annotations:
#'
#' When \code{annotation_row} parameter in \code{\link[pheatmap]{pheatmap}} is used to add row annotations,
#' extra color vectors must be added to the \code{annotation.colors} input to accomodate these.
#' The extra vectors should be named with the annotation bar titles, should come after vectors provided for column annotations.
#' If the vectors are not named, the \code{pheatmap}'s default colors will be used instead.
#' @return A \code{pheatmap} object.
#' @seealso \code{\link[pheatmap]{pheatmap}}, for how to add additional heatmap tweaks.
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = "pbmc",
#'     col.annotation.metas = "ident")
#' #For real data, you will have more cells than this truncated dataset,
#' # so I recommend turning off cell clustering when you are trying out tweaks by
#' # adding cluster_cols=FALSE
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     object = "pbmc",
#'     col.annotation.metas = "ident",
#'     cluster_cols=FALSE)
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' dittoHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'     "FCGR3A","LYZ","PPBP","CD8A"),
#'     col.annotation.metas = "ident")
#'
#' @export

dittoHeatmap <- function(
    genes=NULL, object = DEFAULT, cells.use = NULL, main = NA,
    cell.names.meta = NULL, data.type = "normalized",
    heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
    scaled.to.max = FALSE,
    heatmap.colors.max.scaled = colorRampPalette(c("white", "red"))(25),
    col.annotation.metas = NULL,
    annotation.colors = rep(list(MYcolors),length(col.annotation.metas)),
    data.out=FALSE, highlight.genes = NULL, show.colnames = TRUE,
    show.rownames = TRUE, ...) {

    #If no genes given, "error out"
    if (is.null(genes)) {
        stop('This function is not set up to select which genes to use.\nPlease provide a list of genes a set of genes.')
    }

    #Turn the object into a "name"
    if (typeof(object)=="S4") {
        object <- deparse(substitute(object))
    }

    #If cells.use given as logical, populate as names.
    cells.use <- which_cells(cells.use, object)
    all.cells <- all_cells(object)

    #Make the data matrix
    data <- as.matrix(which_data(data.type,object)[genes,cells.use])
    if (scaled.to.max) {
        maxs <- apply(data,1,max)
        data <- data/maxs
        heatmap.colors <- heatmap.colors.max.scaled
    }

    if (is.null(main)) {
        main <- NA
    }

    args <- list(
        mat = data,
        main = main,
        show_colnames = show.colnames,
        show_rownames = show.rownames,
        color = heatmap.colors,
        scale = ifelse(scaled.to.max, "none", "row"),
        ...)

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
        args$labels_row[-(match(highlight.genes,Row_names))] <- ""
        args$show_colnames <- TRUE
    }

    #Make the columns annotations data for colored annotation bars
    if (!is.null(col.annotation.metas)) {
        args$annotation_col <- data.frame(
            vapply(
                col.annotation.metas,
                function(X) as.character(meta(X,object)[all.cells %in% cells.use]),
                FUN.VALUE = character(length(cells.use))),
            row.names = cells.use)
        if (is.null(names(annotation.colors))) {
            names(annotation.colors)[seq_along(col.annotation.metas)] <- col.annotation.metas
            for (i in seq_along(col.annotation.metas)){
                names(annotation.colors[[i]]) <- levels(as.factor(args$annotation_col[,i]))
                annotation.colors[[i]] <- annotation.colors[[i]][seq_along(levels(as.factor(args$annotation_col[,i])))]
            }
        }
        args$annotation_colors <- annotation.colors
    }

    if (data.out) {
        OUT <- list(args = args, call = "do.call(pheatmap, args)")
    } else {
        OUT <- do.call(pheatmap::pheatmap, args)
    }
    OUT
}
