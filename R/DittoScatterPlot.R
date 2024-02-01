#' Show RNAseq data overlayed on a scatter plot
#' @import ggplot2
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param x.var,y.var Single string giving a gene or metadata that will be used for the x- and y-axis of the scatterplot.
#' Note: must be continuous.
#'
#' Alternatively, can be a directly supplied numeric vector of length equal to the total number of cells/samples in \code{object}.
#' @param color.var Single string giving a gene or metadata that will set the color of cells/samples in the plot.
#'
#' Alternatively, can be a directly supplied numeric or string vector or a factor of length equal to the total number of cells/samples in \code{object}.
#' @param shape.by Single string giving a metadata (Note: must be discrete.) that will set the shape of cells/samples in the plot.
#'
#' Alternatively, can be a directly supplied string vector or a factor of length equal to the total number of cells/samples in \code{object}.
#' @param split.by 1 or 2 strings naming discrete metadata to use for splitting the cells/samples into multiple plots with ggplot faceting.
#'
#' When 2 metadatas are named, c(row,col), the first is used as rows and the second is used for columns of the resulting grid.
#'
#' When 1 metadata is named, shape control can be achieved with \code{split.nrow} and \code{split.ncol}
#'
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}.
#' @param extra.vars String vector providing names of any extra metadata to be stashed in the dataframe supplied to \code{ggplot(data)}.
#'
#' Useful for making custom alterations \emph{after} dittoSeq plot generation.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param show.others Logical. FALSE by default, whether other cells should be shown in the background in light gray.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default, see \code{\link{dittoColors}} for contents.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use.
#'
#' Useful for quickly swapping the colors of nearby clusters.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param assay.x,assay.y,assay.color,assay.extra,slot.x,slot.y,slot.color,slot.extra single strings or integers (SCEs and SEs) or an optionally named vector of such values that set which expression data to use for each given data target.
#' See \code{\link{GeneTargeting}} for specifics and examples -- Seurat and SingleCellExperiment objects deal with these differently, and functionality additions in dittoSeq have led to some minimal divergence from the native methodologies.
#' @param adjustment.x,adjustment.y,adjustment.color,adjustment.extra For the given data target, when targeting gene / feature expression, should that data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param do.hover Logical which controls whether the object will be converted to a plotly object so that data about individual points will be displayed when you hover your cursor over them.
#' \code{hover.data} argument is used to determine what data to use.
#' @param hover.data String vector of gene and metadata names, example: \code{c("meta1","gene1","meta2","gene2")} which determines what data to show on hover when \code{do.hover} is set to \code{TRUE}.
#' @param hover.assay,hover.slot,hover.adjustment Similar to the x, y, color, and extra versions, when showing expression data upon hover, these set what data will be shown.
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.by}, this sets the panel of shapes.
#' When nothing is supplied to \code{shape.by}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#'
#' Note: Unfortunately, shapes can be hard to see when points are on top of each other & they are more slowly processed by the brain.
#' For these reasons, even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param size Number which sets the size of data points. Default = 1.
#' @param opacity Number between 0 and 1.
#' Great for when you have MANY overlapping points, this sets how solid the points should be:
#' 1 = not see-through at all. 0 = invisible. Default = 1.
#' (In terms of typical ggplot variables, = alpha)
#' @param rename.color.groups,rename.shape.groups String vector containing new names for the identities of the color or shape overlay groups.
#' @param add.trajectory.curves List of matrices, each representing coordinates for a trajectory path, from start to end, where matrix columns represent x and y coordinates of the paths.
#' @param legend.show Logical. Whether any legend should be displayed. Default = \code{TRUE}.
#' @param legend.color.title,legend.shape.title Strings which set the title for the color or shape legends.
#' @param legend.color.size,legend.shape.size Numbers representing the size at which shapes should be plotted in the color and shape legends (for discrete variable plotting).
#' Default = 5. *Enlarging the icons in the colors legend is incredibly helpful for making colors more distinguishable by color blind individuals.
#' @param min.color color for \code{min} value of \code{color.var} data. Default = yellow
#' @param max.color color for \code{max} value of \code{color.var} data. Default = blue
#' @param min,max Number which sets the values associated with the minimum or maximum colors.
#' @param legend.color.breaks Numeric vector which sets the discrete values to label in the color-scale legend for continuous data.
#' @param legend.color.breaks.labels String vector, with same length as \code{legend.breaks}, which sets the labels for the tick marks of the color-scale.
#' @param main String, sets the plot title.
#' A default title is automatically generated if based on \code{color.var} and \code{shape.by} when either are provided.
#' To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle.
#' @param xlab,ylab Strings which set the labels for the axes. To remove, set to \code{NULL}.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_bw()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param do.raster Logical. When set to \code{TRUE}, rasterizes the internal plot layer, changing it from individually encoded points to a flattened set of pixels.
#' This can be useful for editing in external programs (e.g. Illustrator) when there are many thousands of data points.
#' @param raster.dpi Number indicating dots/pixels per inch (dpi) to use for rasterization. Default = 300.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot ("p"),
#' a data.frame containing the underlying data for target cells ("Target_data"),
#' and a data.frame containing the underlying data for non-target cells ("Others_data").
#' 
#' @inheritParams dittoDimPlot
#' @return a ggplot scatterplot where colored dots and/or shapes represent individual cells/samples. X and Y axes can be gene expression, numeric metadata, or manually supplied values.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' Alternatively, if \code{do.hover} is set to \code{TRUE}, the plot is coverted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' @details
#' This function creates a dataframe with X, Y, color, shape, and faceting data determined by \code{x.var}, \code{y.var}, \code{color.var}, \code{shape.var}, and \code{split.by}.
#' Any extra gene or metadata requested with \code{extra.var} is added as well.
#' For expression/counts data, \code{assay}, \code{slot}, and \code{adjustment} inputs (\code{.x}, \code{.y}, and \code{.color}) can be used to change which data is used, and if it should be adjusted in some way.
#'
#' Next, if a set of cells or samples to use is indicated with the \code{cells.use} input, then the dataframe is split into \code{Target_data} and \code{Others_data} based on subsetting by the target cells/samples.
#'
#' Finally, a scatter plot is created using these dataframes.
#' Non-target cells are colored in gray if \code{show.others=TRUE},
#' and target cell data is displayed on top, colored and shaped based on the \code{color.var}- and \code{shape.by}-associated data.
#' If \code{split.by} was used, the plot will be split into a matrix of panels based on the associated groupings.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item \code{size} and \code{opacity} can be used to adjust the size and transparency of the data points.
#' \item Colors used can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes used can be adjusted with \code{shape.panel}.
#' \item Color and shape labels can be changed using \code{rename.color.groups} and \code{rename.shape.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab completion lookup.
#' }
#'
#' @seealso
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{x.var}, \code{y.var}, \code{color.var}, \code{shape.by}, and \code{hover.data} options are of an \code{object}.
#'
#' \code{\link{dittoDimPlot}} for making very similar data representations, but where dimensionality reduction (PCA, t-SNE, UMAP, etc.) dimensions are the scatterplot axes.
#'
#' \code{\link{dittoDimHex}} and \code{\link{dittoScatterHex}} for showing very similar data representations, but where nearby cells are summarized together in hexagonal bins.
#'
#' @author Daniel Bunis and Jared Andrews
#' @export
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' # Mock up some nCount_RNA and nFeature_RNA metadata
#' #  == the default way to extract
#' myRNA$nCount_RNA <- runif(60,200,1000)
#' myRNA$nFeature_RNA <- myRNA$nCount_RNA*runif(60,0.95,1.05)
#' # and also percent.mito metadata
#' myRNA$percent.mito <- sample(c(runif(50,0,0.05),runif(10,0.05,0.2)))
#'
#' dittoScatterPlot(
#'     myRNA, x.var = "nCount_RNA", y.var = "nFeature_RNA")
#'
#' # Shapes or colors can be overlaid representing discrete metadata
#' #   or (only colors) continuous metadata / expression data by providing
#' #   metadata or gene names to 'color.var' and 'shape.by'
#' dittoScatterPlot(
#'     myRNA, x.var = "gene1", y.var = "gene2",
#'     color.var = "groups",
#'     shape.by = "SNP",
#'     size = 3)
#' dittoScatterPlot(
#'     myRNA, x.var = "gene1", y.var = "gene2",
#'     color.var = "gene3")
#' 
#' # Note: scatterplots like this can be very useful for dataset QC, especially
#' #   with percentage of mitochondrial reads as the color overlay.
#' dittoScatterPlot(myRNA,
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA",
#'     color.var = "percent.mito")
#'
#' # Data can be "split" or faceted by a discrete variable as well.
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     split.by = "timepoint") # single split.by element
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     split.by = c("groups","SNP")) # row and col split.by elements
#' # OR with 'extra.vars' plus manually faceting for added control
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     extra.vars = c("SNP")) +
#'     facet_wrap("SNP", ncol = 1, strip.position = "left")
#' 
#' # Countours can also be added to help illuminate overlapping samples
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     do.contour = TRUE)
#' 
#' # Multiple continuous metadata or genes can also be plotted together by
#' #   giving that vector to 'color.var':
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     color.var = c("gene3", "gene4"))
#' # This functionality can be combined with 1 additional 'split.by' variable,
#' #   with the directionality then controlled via 'multivar.split.dir':
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     color.var = c("gene3", "gene4"),
#'     split.by = "timepoint",
#'     multivar.split.dir = "col")
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     color.var = c("gene3", "gene4"),
#'     split.by = "timepoint",
#'     multivar.split.dir = "row")
#'
dittoScatterPlot <- function(
    object,
    x.var,
    y.var,
    color.var = NULL,
    shape.by = NULL,
    split.by = NULL,
    extra.vars = NULL,
    cells.use = NULL,
    multivar.split.dir = c("col", "row"),
    show.others = FALSE,
    split.show.all.others = TRUE,
    size = 1,
    opacity = 1,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
    assay.x = .default_assay(object),
    slot.x = .default_slot(object),
    adjustment.x = NULL,
    assay.y = .default_assay(object),
    slot.y = .default_slot(object),
    adjustment.y = NULL,
    assay.color = .default_assay(object),
    slot.color = .default_slot(object),
    adjustment.color = NULL,
    assay.extra = .default_assay(object),
    slot.extra = .default_slot(object),
    adjustment.extra = NULL,
    swap.rownames = NULL,
    shape.panel = c(16,15,17,23,25,8),
    rename.color.groups = NULL,
    rename.shape.groups = NULL,
    min.color = "#F0E442",
    max.color = "#0072B2",
    min = NA,
    max = NA,
    order = c("unordered", "increasing", "decreasing", "randomize"),
    xlab = x.var,
    ylab = y.var,
    main = "make",
    sub = NULL,
    theme = theme_bw(),
    do.hover = FALSE,
    hover.data = NULL,
    hover.assay = .default_assay(object),
    hover.slot = .default_slot(object),
    hover.adjustment = NULL,
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    add.trajectory.lineages = NULL,
    add.trajectory.curves = NULL,
    trajectory.cluster.meta,
    trajectory.arrow.size = 0.15,
    do.letter = FALSE,
    do.ellipse = FALSE,
    do.label = FALSE,
    labels.size = 5,
    labels.highlight = TRUE,
    labels.repel = TRUE,
    labels.split.by = split.by,
    labels.repel.adjust = list(),
    legend.show = TRUE,
    legend.color.title = "make",
    legend.color.size = 5,
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.shape.title = shape.by,
    legend.shape.size = 5,
    do.raster = FALSE,
    raster.dpi = 300,
    data.out = FALSE) {

    order <- match.arg(order)
    multivar.split.dir <- match.arg(multivar.split.dir)
    
    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    # Make dataframe
    all_data <- .scatter_data_gather(
        object = object, cells.use = cells.use, x.var = x.var, y.var = y.var,
        color.var = color.var, shape.by = shape.by, split.by = split.by,
        extra.vars = extra.vars, multivar.split.dir = multivar.split.dir,
        assay.x = assay.x, slot.x = slot.x, adjustment.x = adjustment.x,
        assay.y = assay.y, slot.y = slot.y, adjustment.y = adjustment.y,
        assay.color = assay.color, slot.color = slot.color,
        adjustment.color = adjustment.color,
        assay.extra = assay.extra, slot.extra = slot.extra,
        adjustment.extra = adjustment.extra,
        swap.rownames = swap.rownames,
        do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
        rename.color.groups, rename.shape.groups)
    Target_data <- all_data$Target_data
    Others_data <- all_data$Others_data
    split.by <- all_data$split.by
    
    if (order %in% c("increasing", "decreasing")) {
        Target_data <- Target_data[order(Target_data$color, decreasing = order=="decreasing"),]
    } else if (order == "randomize") {
        Target_data <- Target_data[sample(nrow(Target_data)),]
    }

    # Set title if "make"
    main <- .leave_default_or_null(main,
        paste0(c(color.var, shape.by), collapse = " and "),
        null.if = length(color.var)>1)
    legend.color.title <- .leave_default_or_null(
        legend.color.title, color.var, null.if = length(color.var)>1)

    # Make the plot
    p <- .ditto_scatter_plot(Target_data, Others_data,
        color.var, shape.by, show.others, size, opacity,
        color.panel, colors, do.hover, shape.panel,
        min.color, max.color, min, max,
        xlab, ylab, main, sub, theme,
        legend.show, legend.color.title, legend.color.size,
        legend.color.breaks, legend.color.breaks.labels, legend.shape.title,
        legend.shape.size, do.raster, raster.dpi,
        split.by, split.show.all.others)

    ### Add extra features
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, split.adjust)
    }
    
    if (do.contour) {
        p <- .add_contours(p, Target_data, contour.color, contour.linetype)
    }
    
    p <- .add_letters_ellipses_labels_if_discrete(
        p, Target_data, is.discrete = !is.numeric(Target_data$color),
        do.letter, do.ellipse, do.label,
        labels.highlight, labels.size, labels.repel, labels.split.by,
        labels.repel.adjust,
        size, opacity, legend.color.title, legend.color.size)
    
    if (is.list(add.trajectory.lineages)) {
        p <- .add_trajectory_lineages(
            p, rbind(Target_data, Others_data), add.trajectory.lineages,
            trajectory.cluster.meta, trajectory.arrow.size, object)
    }
    
    if (is.list(add.trajectory.curves)) {
        p <- .add_trajectory_curves(
            p, add.trajectory.curves, trajectory.arrow.size)
    }
    
    if (do.hover) {
        .error_if_no_plotly()
        p <- plotly::ggplotly(p, tooltip = "text")
    }

    ### RETURN the PLOT ###
    if (data.out) {
        list(
            plot = p,
            Target_data = Target_data,
            Others_data = Others_data)
    } else{
        p
    }
}

.ditto_scatter_plot <- function(
    Target_data,
    Others_data,
    color.var,
    shape.by,
    show.others,
    size,
    opacity,
    color.panel,
    colors,
    do.hover,
    shape.panel,
    min.color,
    max.color,
    min,
    max,
    xlab,
    ylab,
    main,
    sub,
    theme,
    legend.show,
    legend.color.title,
    legend.color.size,
    legend.color.breaks,
    legend.color.breaks.labels,
    legend.shape.title,
    legend.shape.size,
    do.raster,
    raster.dpi,
    split.by,
    split.show.all.others
) {
    
    ### Set up plotting
    p <- ggplot() + ylab(ylab) + xlab(xlab) + ggtitle(main,sub) + theme

    # Determine how to add data while adding proper theming
    aes.use <- aes(x = .data$X, y = .data$Y)
    geom.args <- list(
        data = Target_data,
        size=size, alpha = opacity)

    if (!is.null(color.var)) {
        
        aes.use <- modifyList(aes.use, aes(color = .data$color))
        
        if (is.numeric(Target_data$color)) {
            p <- p +
            scale_colour_gradient(
                name = legend.color.title, low= min.color, high = max.color,
                limits = c(min, max),
                breaks = legend.color.breaks,
                labels = legend.color.breaks.labels)
        } else {
            p <- p +
            scale_colour_manual(
                name = legend.color.title,
                values = color.panel[colors]) +
            guides(color = guide_legend(override.aes = list(size=legend.color.size)))
        }
    }

    if (!is.null(shape.by)) {
        
        aes.use <- modifyList(aes.use, aes(shape = .data$shape))
        
        p <- p +
            scale_shape_manual(
                values = shape.panel[seq_along(levels(as.factor(Target_data$shape)))],
                name = legend.shape.title) +
            guides(shape = guide_legend(override.aes = list(size=legend.shape.size)))
    
    } else {
        geom.args$shape <- shape.panel[1]
    }

    ### Add data
    # Others_data
    if (show.others) {
        if (!is.null(split.by) && split.show.all.others) {
            Others_data <- .rep_all_data_per_facet(
                Target_data, Others_data, split.by)
        }
        
        if (nrow(Others_data)>1) {
            if (do.raster) {
                .error_if_no_ggrastr()
                p <- p + ggrastr::geom_point_rast(data = Others_data,
                    aes(x = .data$X, y = .data$Y), size=size, color = "gray90", raster.dpi = raster.dpi)
            } else {
                p <- p + geom_point(data = Others_data,
                    aes(x = .data$X, y = .data$Y), size=size, color = "gray90")
            }
        }
    }
    # Target_data
    if (do.hover) {
        aes.use <- modifyList(aes.use, aes(text = .data$hover.string))
        geom.args$mapping <- aes.use
        p <- p + suppressWarnings(do.call(geom_point, geom.args))
    } else {
        geom.args$mapping <- aes.use
        if (do.raster) {
            .error_if_no_ggrastr()
            p <- p + do.call(ggrastr::geom_point_rast, geom.args)
        } else {
            p <- p + do.call(geom_point, geom.args)
        }
    }

    if (!legend.show) {
        p <- .remove_legend(p)
    }

    p
}

.rep_all_data_per_facet <- function(Target_data, Others_data, split.by) {
    
    all_data <- rbind(Target_data, Others_data)
    
    facet <- if (is.null(split.by)) {
        "filler"
    } else {
        do.call(paste, all_data[,split.by, drop = FALSE])
    }
    
    Others_data <- data.frame(row.names = rownames(all_data))
    
    Others_data <- do.call(
        rbind,
        lapply(
            unique(facet),
            function(this_facet) {
        
                facet_data <- all_data[facet==this_facet, , drop = FALSE]
                
                new_data <- all_data        
                # Add facet info
                if (!is.null(split.by)) {
                    for (by in split.by) {
                        new_data[[by]] <- facet_data[1,by]
                    }
                }
                
                new_data
            }
        )
    )
}

.scatter_data_gather <- function(
    object,
    cells.use,
    x.var,
    y.var,
    color.var,
    shape.by,
    split.by,
    extra.vars,
    multivar.split.dir,
    assay.x,
    slot.x,
    adjustment.x,
    assay.y,
    slot.y,
    adjustment.y,
    assay.color,
    slot.color,
    adjustment.color,
    assay.extra,
    slot.extra,
    adjustment.extra,
    swap.rownames = NULL,
    do.hover = FALSE,
    hover.data = NULL,
    hover.assay = NULL,
    hover.slot = NULL,
    hover.adjustment = NULL,
    rename.color.groups = NULL,
    rename.shape.groups = NULL
    ) {

    all.cells <- .all_cells(object)
    cells.use <- .which_cells(cells.use, object)
    object <- .swap_rownames(object, swap.rownames)
    
    # Support multiple genes/metadata
    if (length(color.var)>1 && length(color.var) != length(all.cells)) {
        # Data
        each_data <- lapply(
            color.var, function(this.var) {
                this.out <- .scatter_data_gather_inner(
                    object, all.cells, x.var, y.var, this.var,
                    shape.by, split.by, extra.vars, assay.x, slot.x, adjustment.x,
                    assay.y, slot.y, adjustment.y, assay.color, slot.color,
                    adjustment.color, assay.extra, slot.extra, adjustment.extra,
                    do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
                    rename.color.groups, rename.shape.groups)
                this.out$color.var <- this.var
                this.out
            }
        )
        if (any(unlist(lapply(each_data, function(x) { !is.numeric(x$color) })))) {
            stop("Only numeric data supported when plotting multiple color.var")
        }
        out_data <-list(
            Target_data = do.call(rbind, lapply(each_data, function(x) { x[cells.use,] } )),
            Others_data = do.call(rbind, lapply(each_data, function(x) { x[!(all.cells %in% cells.use),] } ))
        )
        split.by <- .multivar_adjust_split_by(
            split.by, multivar.split.dir, multivar.col.name = "color.var")
    } else {
        # Single var
        all_data <- .scatter_data_gather_inner(
            object, all.cells, x.var, y.var, color.var,
            shape.by, split.by, extra.vars, assay.x, slot.x, adjustment.x,
            assay.y, slot.y, adjustment.y, assay.color, slot.color,
            adjustment.color, assay.extra, slot.extra, adjustment.extra,
            do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
            rename.color.groups, rename.shape.groups)
        out_data <-list(
            Target_data = all_data[cells.use,],
            Others_data = all_data[!(all.cells %in% cells.use),]
        )
    }
    
    out_data$split.by <- split.by
    
    out_data
}

.scatter_data_gather_inner <- function(
    object,
    all.cells,
    x.var,
    y.var,
    color.var,
    shape.by,
    split.by,
    extra.vars,
    assay.x,
    slot.x,
    adjustment.x,
    assay.y,
    slot.y,
    adjustment.y,
    assay.color,
    slot.color,
    adjustment.color,
    assay.extra,
    slot.extra,
    adjustment.extra,
    do.hover,
    hover.data,
    hover.assay,
    hover.slot,
    hover.adjustment,
    rename.color.groups,
    rename.shape.groups) {
    
    # Make dataframe
    vars <- list(x.var, y.var, color.var, shape.by)
    names <- list("X", "Y", "color", "shape")
    assays <- list(assay.x, assay.y, assay.color, NA)
    slots <- list(slot.x, slot.y, slot.color, NA)
    adjustments <- list(adjustment.x, adjustment.y, adjustment.color, NA)
    relabels <- list(NULL, NULL, rename.color.groups, rename.shape.groups)

    dat <- data.frame(row.names = all.cells)
    for (i in seq_along(vars)) {
        dat <- .add_by_cell(dat, vars[[i]], names[[i]], object, assays[[i]],
            slots[[i]], adjustments[[i]], NULL, relabels[[i]])
    }

    extra.vars <- unique(c(split.by, extra.vars))
    dat <- .add_by_cell(dat, extra.vars, extra.vars, object, assay.extra,
        slot.extra, adjustment.extra, mult = TRUE)

    if (do.hover) {
        dat$hover.string <- .make_hover_strings_from_vars(
            hover.data, object, hover.assay, hover.slot, hover.adjustment)
    }
    
    dat
}

.multivar_adjust_split_by <- function(split.by, multivar.split.dir, multivar.col.name) {
    if (is.null(split.by)) {
        split.by <- multivar.col.name
    } else {
        if (length(split.by)>1) {
            warning("multi-feature '", multivar.col.name, "' is prioiritized for faceting. The second 'split.by' element will be ignored.")
        }
        split.by[2] <- multivar.col.name
        if (multivar.split.dir=="row") {
            split.by <- rev(split.by)
        }
    }
    split.by
}
