#' Shows Neighbor Diversity data, per a given metadata, overlaid on a umap, tsne, pca, or similar
#' @name dittoNeighborDiversity
#' @import ggplot2
#'
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample. Discrete or continuous data both work.
#'
#' Alternatively, a string vector naming multiple genes or metadata, OR a vector of the same length as there are cells/samples in the \code{object} which provides per-cell data directly.
#' @param size Number which sets the size of data points. Default = 0.1 here to enable seeing more cells in dense regions.
#' @param opacity Number between 0 and 1, which defaults to 0.8 here, and can be increased or lowered to make cells less or more transparent, respectively.
#' @param min Number which sets the value associated with the minimum color.  Defaults to 0 here.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing
#' the calculated neighborhood diversity metadata ("diversity") either as vector or data.frame depending on how many metadata were given to \code{var},
#' the plot ("p"),
#' a data.frame containing the underlying data for target cells ("Target_data"),
#' and a data.frame containing the underlying data for non-target cells ("Others_data").
#' @inheritParams dittoDimPlot
#' @inheritParams dittoHex
#'
#' @return A ggplot or plotly object where neighborhood diversity of \code{var}-values among cells' 'nearby' nearest neighbors is overlaid, via color, onto a tSNE, PCA, UMAP, ..., plot of choice.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing four slots is output:
#' the calculated neighborhood diversity metadata ("diversity") either as vector or data.frame depending on how many metadata were given to \code{var},
#' the plot (named 'p'),
#' a data.table containing the underlying data for target cells (named 'Target_data'),
#' and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' Alternatively, if \code{do.hover} is set to \code{TRUE}, the plot is coverted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' @details
#' These plotters start by making use of \code{\link{calcNeighborMetadataDiversity}}, then
#' passes all inputs through to \code{\link{dittoDimPlot}} or \code{\link{dittoDimHex}} plotters
#' 
#' @seealso
#' \code{\link{calcNeighborMetadataDiversity}} for details on the neighborhood diversity calculation
#' \code{\link{dittoDimPlot}} and \code{\link{dittoDimHex}} for additional details about the other options as these are the plotters used after diversity calculations complete.
#'
#' @author Daniel Bunis
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#' 
#' # Temporary Seurat code for calculating neighbors
#' dittoSeq:::.error_if_no_Seurat()
#' myRNA <- Seurat::as.Seurat(myRNA)
#' myRNA <- Seurat::FindNeighbors(myRNA, reduction = "pca", dims = 1:5, return.neighbor = TRUE)
#' 
#' # (Using bigger size than the default for these examples because the example data has so few cells)
#' dittoNeighborDiversityPlot(myRNA, "groups", size = 1)
#' dittoNeighborDiversityHex(myRNA, "groups")
#'
NULL

#' @describeIn dittoNeighborDiversity Shows Neighbor Diversity data, per a given metadata, overlaid per cell on a umap, tsne, pca, or similar
#' @export
dittoNeighborDiversityPlot <- function(
    object,
    var,
    neighbors = .default_neighbors(object),
    distances,
    quantile = 0.9,
    reduction.use = .default_reduction(object),
    size = 0.1,
    opacity = 1,
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    shape.by = NULL,
    split.by = NULL,
    split.adjust = list(),
    extra.vars = NULL,
    multivar.split.dir = c("col", "row"),
    show.others = TRUE,
    split.show.all.others = TRUE,
    split.nrow = NULL,
    split.ncol = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    shape.panel = c(16,15,17,23,25,8),
    min.color = "#F0E442",
    max.color = "#0072B2",
    min = NA,
    max = NA,
    order = c("unordered", "increasing", "decreasing", "randomize"),
    main = "make",
    sub = NULL,
    xlab = "make",
    ylab = "make",
    rename.var.groups = NULL,
    rename.shape.groups = NULL,
    theme = theme_bw(),
    show.axes.numbers = TRUE,
    show.grid.lines = if (is.character(reduction.use)) { !grepl("umap|tsne", tolower(reduction.use)) } else {TRUE},
    do.hover = FALSE,
    hover.data = c(var, paste0(var, "_diversity")),
    add.trajectory.lineages = NULL,
    add.trajectory.curves = NULL,
    trajectory.cluster.meta,
    trajectory.arrow.size = 0.15,
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    legend.show = TRUE,
    legend.size = 5,
    legend.title = "make",
    legend.breaks = waiver(),
    legend.breaks.labels = waiver(),
    shape.legend.size = 5,
    shape.legend.title = shape.by,
    do.raster = FALSE,
    raster.dpi = 300,
    data.out = FALSE) {
    
    var_use <- c()
    for (this_var in var) {
        this_var_use <- paste0(this_var, "_diversity")
        object[[this_var_use]] <- calcNeighborMetadataDiversity(
            object, this_var, neighbors, distances, quantile
        )
        var_use <- c(var_use, this_var_use)
    }
    div_out <- getMetas(object, names.only = FALSE)[,var_use]
    
    # Make dataframes and plot
    p.df <- dittoDimPlot(
        object, var_use, reduction.use, size, opacity, dim.1, dim.2, cells.use,
        shape.by, split.by, split.adjust, extra.vars, multivar.split.dir,
        show.others, split.show.all.others, split.nrow, split.ncol,
        assay = NA, slot = NA, adjustment = NULL, swap.rownames = NULL,
        color.panel, colors, shape.panel, min.color, max.color, min, max, order,
        main, sub, xlab, ylab, rename.var.groups, rename.shape.groups, theme,
        show.axes.numbers, show.grid.lines,
        do.letter = FALSE, do.ellipse = FALSE, do.label = FALSE,
        labels.size = 5, labels.highlight = TRUE, labels.repel = TRUE,
        labels.split.by = split.by, labels.repel.adjust = list(),
        do.hover, hover.data,
        hover.assay = NA, hover.slot = NA, hover.adjustment = NULL,
        add.trajectory.lineages, add.trajectory.curves, trajectory.cluster.meta,
        trajectory.arrow.size, do.contour, contour.color, contour.linetype,
        legend.show, legend.size, legend.title, legend.breaks,
        legend.breaks.labels, shape.legend.size, shape.legend.title,
        do.raster, raster.dpi, data.out = TRUE)
    p <- p.df$plot
    Target_data <- p.df$Target_data
    Others_data <- p.df$Others_data
    
    ### RETURN the PLOT ###
    if (data.out) {
        list(
            diversity = div_out,
            plot = p,
            Target_data = Target_data,
            Others_data = Others_data)
    } else {
        p
    }
}

#' @describeIn dittoNeighborDiversity Shows Neighbor Diversity data, per a given metadata, summarized and overlaid per hexagonaly-shaped region on a umap, tsne, pca, or similar plot
#' @export
dittoNeighborDiversityHex <- function(
        object,
        var,
        neighbors = .default_neighbors(object),
        distances,
        quantile = 0.9,
        bins = 30,
        color.method = NULL,
        reduction.use = .default_reduction(object),
        dim.1 = 1,
        dim.2 = 2,
        cells.use = NULL,
        color.panel = dittoColors(),
        colors = seq_along(color.panel),
        split.by = NULL,
        extra.vars = NULL,
        multivar.split.dir = c("col", "row"),
        split.nrow = NULL,
        split.ncol = NULL,
        split.adjust = list(),
        assay.extra = assay,
        slot.extra = slot,
        adjustment.extra = adjustment,
        show.axes.numbers = TRUE,
        show.grid.lines = !grepl("umap|tsne", tolower(reduction.use)),
        main = "make",
        sub = NULL,
        xlab = "make",
        ylab = "make",
        theme = theme_bw(),
        do.contour = FALSE,
        contour.color = "black",
        contour.linetype = 1,
        min.density = NA,
        max.density = NA,
        min.color = "#F0E442",
        max.color = "#0072B2",
        min.opacity = 0.2,
        max.opacity = 1, 
        min = NA,
        max = NA,
        do.ellipse = FALSE,
        add.trajectory.lineages = NULL,
        add.trajectory.curves = NULL,
        trajectory.cluster.meta,
        trajectory.arrow.size = 0.15,
        data.out = FALSE,
        legend.show = TRUE,
        legend.color.title = "make",
        legend.color.breaks = waiver(),
        legend.color.breaks.labels = waiver(),
        legend.density.title = if (isBulk(object)) "Samples" else "Cells",
        legend.density.breaks = waiver(),
        legend.density.breaks.labels = waiver()
) {
    var_use <- c()
    for (this_var in var) {
        this_var_use <- paste0(this_var, "_diversity")
        object[[this_var_use]] <- calcNeighborMetadataDiversity(
            object, this_var, neighbors, distances, quantile
        )
        var_use <- c(var_use, this_var_use)
    }
    div_out <- getMetas(object, names.only = FALSE)[,var_use]
    
    # Make dataframes and plot
    p.df <- dittoDimHex(
        object, var_use, bins, color.method, reduction.use, dim.1, dim.2,
        cells.use, color.panel, colors, split.by, extra.vars,
        multivar.split.dir, split.nrow, split.ncol, split.adjust,
        assay = NA, slot = NA, adjustment = NULL, swap.rownames = NULL,
        assay.extra = assay, slot.extra = slot, adjustment.extra = adjustment,
        show.axes.numbers, show.grid.lines, main, sub, xlab, ylab, theme,
        do.contour, contour.color, contour.linetype, min.density, max.density,
        min.color, max.color, min.opacity, max.opacity, min, max,
        rename.color.groups = NULL, do.ellipse, do.label = FALSE,
        labels.size = 5, labels.highlight = TRUE, labels.repel = TRUE,
        labels.split.by = split.by, labels.repel.adjust = list(),
        add.trajectory.lineages, add.trajectory.curves, trajectory.cluster.meta,
        trajectory.arrow.size, data.out = TRUE, legend.show, legend.color.title,
        legend.color.breaks, legend.color.breaks.labels, legend.density.title,
        legend.density.breaks, legend.density.breaks.labels
        )
    p <- p.df$plot
    data <- p.df$data
    
    ### RETURN the PLOT ###
    if (data.out) {
        list(
            diversity = div_out,
            plot = p,
            data = data)
    } else {
        p
    }
}

#' Calculates Neighbor Diversity of a given metadata per each cell 
#'
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample. Discrete or continuous data both work.
#'
#' Alternatively, a string vector naming multiple genes or metadata, OR a vector of the same length as there are cells/samples in the \code{object} which provides per-cell data directly.
#' @param neighbors a single string giving either the name of a Neighbors slot of the (Seurat) \code{object},
#' OR or matrix with cells in its rows and indexes of neighbors in its columns
#' @param distances not needed when \code{neighbors} is directed to a Neighbors slot of the (Seurat) \code{object},
#' Otherwise, must be given a matrix with cells in its rows and distance measures to each neighbor in its columns
#' @return A named numeric vector of diversity counts, the same length as the number of cells in \code{object} which can be added to the object as cell metadata, and named by the cell names of the \code{object}.
#'
#' @details
#' If given a Seurat \code{object} and \code{neighbors} is given (default) a string value representing a Neighbors object slot name.
#' It then extracts the \code{neighbors}-matrix and \code{distances}-matrix from this object.
#' 
#' Otherwise, it uses the \code{neighbors} and \code{distances} inputs for these purposes.
#' 
#' To calculate neighbors' Diversity:
#' 
#' First, the distance cutoff for neighbors deemed close-enough is determined based on the given \code{quantile} of \code{distance}-values.
#' 
#' Then, function then extracts the given \code{var} metadata from the object.
#' 
#' Finally, it loops through each cell (row) of the neighbors and distances matrices,
#' totaling the number of distinct var-values associated with the cell's neighbors that are within the threshold distance.  
#' 
#' Cell names are then added directly before the vector is output.
#' 
#' @seealso
#' \code{\link{dittoDimPlot}} and \code{\link{dittoDimHex}} for additional details about the other options as these are the plotters used after diversity calculations complete.
#'
#' @author Daniel Bunis
#' @export
#' @examples
#' example(importDittoBulk, echo = FALSE)
calcNeighborMetadataDiversity <- function(object, var, neighbors = .default_neighbors(object), distances, quantile = 0.9) {
    
    if (is.character(neighbors)) {
        .error_if_no_Seurat()
        neighbor_object <- SeuratObject::Neighbors(object, neighbors)
        neighbors <- neighbor_object@nn.idx
        distances <- neighbor_object@nn.dist
    }
    if (!nrow(neighbors)==ncol(object)) {
        stop("The number of cells in 'object' does not match the number of cells tracked in the given 'neighbors' data.")
    }
    if (!nrow(distances)==ncol(object)) {
        stop("The number of cells in 'object' does not match the number of cells tracked in the given 'distances' data.")
        
    }
    threshold <- quantile(distances, 0.9)
    
    compar <- meta(var, object)
    
    OUT <- vapply(
        seq_len(ncol(object)),
        function(i) {
            length(unique(compar[neighbors[i,distances[i,]<=threshold]]))
        },
        double(1))
    names(OUT) <- .all_cells(object)
    OUT
}
