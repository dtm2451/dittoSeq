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
#' @param assay.x,assay.y,assay.color,assay.extra,slot.x,slot.y,slot.color,slot.extra,adjustment.x,adjustment.y,adjustment.color,adjustment.extra assay, slot, and adjustment set which data to use when the axes, coloring, or \code{extra.vars} are based on expression/counts data. See \code{\link{gene}} for additional information.
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
#' Note: \code{do.hover} plotly conversion is turned off in this setting, but hover.data is still calculated.
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
#' # Countours can also be added to help illumunate overlapping samples
#' dittoScatterPlot(myRNA, x.var = "gene1", y.var = "gene2",
#'     do.contour = TRUE)
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
    show.others = FALSE,
    size = 1,
    opacity = 1,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    split.nrow = NULL,
    split.ncol = NULL,
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
    min = NULL,
    max = NULL,
    order = c("unordered", "increasing", "decreasing"),
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
    legend.show = TRUE,
    legend.color.title = color.var,
    legend.color.size = 5,
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.shape.title = shape.by,
    legend.shape.size = 5,
    do.raster = FALSE,
    raster.dpi = 300,
    data.out = FALSE) {

    order <- match.arg(order)
    
    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    # Make dataframe
    all_data <- .scatter_data_gather(
        object = object, cells.use = cells.use, x.var = x.var, y.var = y.var,
        color.var = color.var, shape.by = shape.by, split.by = split.by,
        extra.vars = extra.vars,
        assay.x = assay.x, slot.x = slot.x, adjustment.x = adjustment.x,
        assay.y = assay.y, slot.y = slot.y, adjustment.y = adjustment.y,
        assay.color = assay.color, slot.color = slot.color,
        adjustment.color = adjustment.color,
        assay.extra = assay.extra, slot.extra = slot.extra,
        adjustment.extra = adjustment.extra,
        swap.rownames = swap.rownames,
        do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
        rename.color.groups, rename.shape.groups)

    # Trim by cells.use, then order if wants
    Target_data <- all_data[cells.use,]
    Others_data <- all_data[!(all.cells %in% cells.use),]
    
    if (order != "unordered") {
        decreasing <- switch(order,
            "decreasing" = TRUE,
            "increasing" = FALSE)
        Target_data <- Target_data[order(Target_data$color, decreasing = decreasing),]
    }

    # Set title if "make"
    main <- .leave_default_or_null(main,
        paste0(c(color.var, shape.by), collapse = " and "))

    # Make the plot
    p <- .ditto_scatter_plot(Target_data, Others_data,
        color.var, shape.by, show.others, size, opacity,
        color.panel, colors, do.hover, shape.panel,
        min.color, max.color, min, max,
        xlab, ylab, main, sub, theme,
        legend.show, legend.color.title, legend.color.size,
        legend.color.breaks, legend.color.breaks.labels, legend.shape.title,
        legend.shape.size, do.raster, raster.dpi)

    ### Add extra features
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, object, cells.use)
    }
    
    if (do.contour) {
        p <- .add_contours(p, Target_data, contour.color, contour.linetype)
    }
    
    p <- .add_letters_ellipses_labels_if_discrete(
        p, Target_data, is.discrete = !is.numeric(Target_data$color),
        do.letter, do.ellipse, do.label,
        labels.highlight, labels.size, labels.repel, labels.split.by,
        size, opacity, legend.color.title, legend.color.size)
    
    if (is.list(add.trajectory.lineages)) {
        p <- .add_trajectory_lineages(
            p, all_data, add.trajectory.lineages,
            trajectory.cluster.meta, trajectory.arrow.size, object)
    }
    
    if (is.list(add.trajectory.curves)) {
        p <- .add_trajectory_curves(
            p, add.trajectory.curves, trajectory.arrow.size)
    }

    ### RETURN the PLOT ###
    if (data.out) {
        return(list(plot = p, Target_data = Target_data, Others_data = Others_data))
    } else{
        if (do.hover) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
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
    raster.dpi
) {
    
    ### Set up plotting
    p <- ggplot() + ylab(ylab) + xlab(xlab) + ggtitle(main,sub) + theme

    # Determine how to add data while adding proper theming
    aes.args <- list(x = "X", y = "Y")
    geom.args <- list(
        data = Target_data,
        size=size, alpha = opacity)

    if (!is.null(color.var)) {
        
        aes.args$color = "color"
        
        if (is.numeric(Target_data$color)) {
            p <- p +
            scale_colour_gradient(
                name = legend.color.title, low= min.color, high = max.color,
                limits = c(
                    ifelse(is.null(min), min(Target_data$color), min),
                    ifelse(is.null(max), max(Target_data$color), max)),
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
        
        aes.args$shape = "shape"
        
        p <- p +
            scale_shape_manual(
                values = shape.panel[seq_along(levels(Target_data$shape))],
                name = legend.shape.title) +
            guides(shape = guide_legend(override.aes = list(size=legend.shape.size)))
    
    } else {
        geom.args$shape <- shape.panel[1]
    }

    ### Add data
    if (show.others && nrow(Others_data)>1) {
        if (do.raster) {
            .error_if_no_ggrastr()
            p <- p + ggrastr::geom_point_rast(data = Others_data,
                aes_string(x = "X", y = "Y"), size=size, color = "gray90", raster.dpi = raster.dpi)
        } else {
            p <- p + geom_point(data = Others_data,
                aes_string(x = "X", y = "Y"), size=size, color = "gray90")
        }
    }

    if (do.hover) {
        aes.args$text = "hover.string"
        geom.args$mapping <- do.call(aes_string, aes.args)
        p <- p + suppressWarnings(do.call(geom_point, geom.args))
    } else {
        geom.args$mapping <- do.call(aes_string, aes.args)
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
