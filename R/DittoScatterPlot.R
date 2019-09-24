#' Show RNAseq data overlayed on a scatter plot
#' @import ggplot2
#'
#' @param x.var Variable for setting x-axis position of cells/samples.  Can be the name of a gene, meta-data, or "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param y.var Variable for setting y-axis position of cells/samples.  Can be the name of a gene, meta-data, or "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param color.var Variable for setting the color of cells/samples in the plot.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param shape.var Variable for setting the shape of cells/samples in the plot.  Note: must be discrete.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param show.others TRUE/FALSE. TRUE by default, whether other cells should be shown in the background
#' @param color.panel a list of colors to be used for when plotting a discrete var.
#' @param colors indexes / order of colors from color.panel to use. USAGE= changing the order of how colors are linked to specific groups
#' @param data.type.x For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param data.type.y For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param data.type.color For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param do.hover TRUE/FALSE. Default = FALSE.  If set to true: object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'hover.data' argument is used to determine what data to use.  NOTE: incompatible with lettering (due to a ggplotly incompatibility). Setting do.hover to TRUE will override a do.letter=TRUE or NA.
#' @param hover.data list of variable names, c("meta1","gene1","meta2","gene2"). determines what data to show on hover when do.hover is set to TRUE.
#' @param hover.data.type For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized" For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.var}, this sets the panel of shapes.
#' When nothing is supplied to \code{shape.var}, only the first value is used.
#' Default is a list of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#' Unfortunately, shapes can be hard to see, especially when points are on top of each other, and they are more slowly processed by the brain.
#' For these reasons, even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param size Number. Size of data points.  Default = 1.
#' @param opacity Number between 0 and 1. Great for when you have MANY overlapping points, this sets how see-through the points should be; 1 = not at all; 0 = invisible. Default = 1.
#' @param rename.color.groups Character vector containing new names for the identities of the color overlay.
#' @param rename.shape.groups Character vector containing new names for the identities of the shape overlay.
#' @param legend.show TRUE/FALSE. Whether the legend should be displayed. Default = TRUE.
#' @param legend.color.title String which sets the title for the color legend. Default = \code{NULL} OR \code{var} when a shape legend will also be shown.
#' @param legend.color.size Number representing the size to increase the plotting of color legend shapes to (for discrete variable plotting).
#' Default = 5. *Enlarging the colors legend is incredibly helpful for making colors more distinguishable by color blind individuals.*
#' @param legend.shape.title For adding a title for the colors legend.  Set to \code{NULL} to turn off
#' @param legend.shape.size Override for legend size of the shape variable.
#' @param min.color color for lowest values shown of the color overlay.  Default is \code{dittoColors()[4]}, a yellow.
#' @param max.color color for highest values shown of the color overlay.  Default is \code{dittoColors()[5]}, a dark blue.
#' @param min set the value associated with the minimum color.  All points with a lower value than this will get the same min.color.
#' @param max set the value associated with the maximum color.  All points with a higher value than this will get the same max.color.  Note: if your legend is not plotting, it may be because min > max.
#' @param legend.color.breaks Numeric vector. Sets the discrete values to show in the color-scale legend for continuous data.
#' @param legend.color.breaks.labels String vector with same length as \code{breaks}. Renames the values displayed next to the color-scale.
#' @param main String, sets the plot title. Default = "make" and if left as make, a title will be automatically generated.  To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle
#' @param xlab String, label for y axes. Default = \code{x.var}. To remove, set to NULL.
#' @param ylab String, label for y axes. Default = \code{y.var}. To remove, set to NULL.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments. Default = \code{theme_bw()}. See \code{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options.
#' @param data.out Whether just the plot should be output, or a list with the plot and Target_data and Others_data dataframes.  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @return Makes a plot where colored dots and/or shapes representing individual cells/samples are overlayed onto a scatterplot where x and y can be gene expression (or any numeric metadata) of those cells/samples.
#' @export
#' @examples
#'
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' # Mock up some percent.mito metadata
#' pbmc$percent.mito <- sample(c(runif(75,0,0.05),runif(5,0.05,0.2)))
#'
#' dittoScatterPlot(
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA", object = "pbmc")
#'
#' # Shapes or colors can be overlaid representing discrete metadata
#' #   or (only colors) continuous metadata / expression data by providing
#' #   metadata or gene names to 'color.var' and 'shape.var'
#' dittoScatterPlot(
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA", object = "pbmc",
#'     color.var = "RNA_snn_res.1",
#'     shape.var = "RNA_snn_res.0.8")
#' dittoScatterPlot(
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA", object = "pbmc",
#'     color.var = "percent.mito",
#'     shape.var = "RNA_snn_res.0.8")
#' dittoScatterPlot(
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA", object = "pbmc",
#'     color.var = "CD14",
#'     shape.var = "RNA_snn_res.0.8")
#'
#' # Note: scatterplots like this can be very useful for dataset QC, epecially
#' #   with percentage of reads coming from genes as the color overlay.
dittoScatterPlot <- function(
    x.var, y.var, color.var = NULL, shape.var = NULL,
    object = DEFAULT, cells.use = NULL, show.others = FALSE,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    size = 1, opacity = 1,
    data.type.x = "normalized", data.type.y = "normalized",
    data.type.color = "normalized", do.hover = FALSE, hover.data = NULL,
    hover.data.type = "normalized", shape.panel=c(16,15,17,23,25,8),
    rename.color.groups = NULL, rename.shape.groups = NULL,
    min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
    xlab = x.var, ylab = y.var, main = NULL, sub = NULL, theme = theme_bw(),
    legend.show = TRUE,
    legend.color.title = color.var, legend.color.size = 5,
    legend.color.breaks = waiver(), legend.color.breaks.labels = waiver(),
    legend.shape.title = shape.var, legend.shape.size = 5,
    data.out = FALSE) {
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    # Make dataframe
    dat <- data.frame(
        X = .var_OR_get_meta_or_gene(x.var, object, data.type.x),
        Y = .var_OR_get_meta_or_gene(y.var, object, data.type.y),
        row.names = all.cells)
    aes.args <- list(x = "X", y = "Y")
    do.color <- FALSE
    do.shape <- FALSE
    if (!is.null(color.var)) {
        dat$color <- .var_OR_get_meta_or_gene(
            color.var, object, data.type.color)
        dat$color <- .rename_and_or_reorder(
            dat$color, relabels = rename.color.groups)
        aes.args$color = "color"
        do.color <- TRUE
    }
    if (!is.null(shape.var)) {
        dat$shape <- .var_OR_get_meta_or_gene(
            shape.var, object, data.type.color)
        dat$shape <- .rename_and_or_reorder(
            dat$shape, relabels = rename.shape.groups)
        aes.args$shape = "shape"
        do.shape <- TRUE
    }
    if (do.hover) {
        hover.string <- .make_hover_strings_from_vars(hover.data, object, hover.data.type)
        aes.args$text = "hover.string"
    }
    Target_data <- dat[cells.use,]
    Others_data <- dat[!(all.cells %in% cells.use),]

    ###Start building the plot###
    p <- ggplot() + ylab(ylab) + xlab(xlab) + ggtitle(main,sub) + theme
    if (do.shape) {
        p <- p +
        scale_shape_manual(
            values = shape.panel[seq_along(levels(Target_data$shape))],
            name = legend.shape.title) +
        guides(shape = guide_legend(override.aes = list(size=legend.shape.size)))
    }
    if (do.color) {
        if (is.numeric(dat$color)) {
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

    ### Add data ###
    if (show.others && nrow(Others_data)>1) {
        p <- p + geom_point(data = Others_data,
            aes_string(x = "X", y = "Y"), size=size, color = "gray90")
    }
    geom.args <- list(
        data = Target_data,
        mapping = do.call(aes_string, aes.args),
        size=size, alpha = opacity)
    if (!do.shape) {
        geom.args$shape <- shape.panel[1]
    }
    p <- p + do.call(geom_point, geom.args)
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    ### RETURN the PLOT ###
    if (data.out) {
        return(list(plot = p, Target_data = Target_data, Others_data = Others_data))
    } else{
        if (do.hover) {
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}
