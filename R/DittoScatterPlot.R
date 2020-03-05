#' Show RNAseq data overlayed on a scatter plot
#' @import ggplot2
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param x.var Single string or numeric vector which sets x-axis position of cells/samples.
#' Note: must be continuous.
#' Can be the name of a gene or metadata.
#' Alternatively, can be a numeric vector of length equal to the total number of cells/samples in object.
#' @param y.var Single string or numeric vector which sets y-axis position of cells/samples.
#' Note: must be continuous.
#' Can be the name of a gene or metadata.
#' Alternatively, can be a numeric vector of length equal to the total number of cells/samples in object.
#' @param color.var Single string or numeric vector which sets the color of cells/samples in the plot.
#' Can be the name of a gene, metadata, or "ident" for clusters of a Seurat object.
#' Alternatively, can be a numeric or string vector or a factor of length equal to the total number of cells/samples in object.
#' @param shape.var Single string or numeric vector which sets the shape of cells/samples in the plot.
#' Note: must be discrete.
#' Can be the name of a metadata, or "ident" for clusters of a Seurat object.
#' Alternatively, can be a string vector or a factor of length equal to the total number of cells/samples in object.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param show.others Logical. TRUE by default, whether other cells should be shown in the background in light gray.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default, see \code{\link{dittoColors}} for contents.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param assay.x,assay.y,assay.color,slot.x,slot.y,slot.color,adjustment.x,adjustment.y,adjustment.color assay, slot, and adjustment set which data to use when the axes or coloring are based on expression data. See \code{\link{gene}} for additional information.
#' @param do.hover Logical which controls whether the object will be converted to a plotly object so that data about individual points will be displayed when you hover your cursor over them.
#' \code{hover.data} argument is used to determine what data to use.
#' @param hover.data String vector of gene and metadata names, example: \code{c("meta1","gene1","meta2","gene2")} which determines what data to show on hover when \code{do.hover} is set to \code{TRUE}.
#' @param hover.assay,hover.slot,hover.adjustment Similar to the x, y, and color versions, when showing expression data upon hover, these set what data will be shown.
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.var}, this sets the panel of shapes.
#' When nothing is supplied to \code{shape.var}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#'
#' Note: Unfortunately, shapes can be hard to see when points are on top of each other & they are more slowly processed by the brain.
#' For these reasons, even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param size Number which sets the size of data points. Default = 1.
#' @param opacity Number between 0 and 1.
#' Great for when you have MANY overlapping points, this sets how solid the points should be:
#' 1 = not see-through at all. 0 = invisible. Default = 1.
#' (In terms of typical ggplot variables, = alpha)
#' @param rename.color.groups String vector containing new names for the identities of the color overlay groups.
#' @param rename.shape.groups String vector containing new names for the identities of the shape overlay groups.
#' @param legend.show Logical. Whether any legend should be displayed. Default = \code{TRUE}.
#' @param legend.color.title,legend.shape.title Strings which sets the title for the color or shape legends. Default = \code{NULL} OR \code{color.var} when a shape legend will also be shown.
#' @param legend.color.size,legend.shape.size Numbers representing the size to increase the plotting of color or shape legend icons to (for discrete variable plotting).
#' Default = 5. *Enlarging the icons in the colors legend is incredibly helpful for making colors more distinguishable by color blind individuals.
#' @param min.color color for lowest values of var/min.  Default = yellow
#' @param max.color color for highest values of var/max.  Default = blue
#' @param min Number which sets the value associated with the minimum color.
#' @param max Number which sets the value associated with the maximum color.
#' @param legend.color.breaks Numeric vector which sets the discrete values to show in the color-scale legend for continuous data.
#' @param legend.color.breaks.labels String vector, with same length as \code{legend.breaks}, which renames what's displayed next to the tick marks of the color-scale.
#' @param main String, sets the plot title.
#' A default title is automatically generated if based on \code{color.var} and \code{shape.var} when either are provided.
#' To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle.
#' @param xlab,ylab Strings which set the labels for the axes. To remove, set to \code{NULL}.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments. Default = \code{theme_bw()}. See \code{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options.
#' @param data.out Whether just the plot should be output, or a list with the plot and Target_data and Others_data dataframes.  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @return Makes a scatterplot from (sc)RNAseq data where colored dots and/or shapes represent individual cells/samples.  X and Y can be gene expression (or any numeric metadata) of those cells/samples.
#' @details
#' This function creates a dataframe with the X and Y coordinates determined by \code{x.var} and \code{y.var}.
#' It then adds data for how coloring should be set if a \code{color.var} is given & for how shapes should be set if a \code{shape.var} is given.
#' The \code{assay}, \code{slot}, and \code{adjustment} inputs (\code{.x}, \code{.y}, and \code{.color}) can be used to change what expression data is used when the target data is gene expression data.
#'
#' Next, if a set of cells or samples to use is indicated with the \code{cells.use} input, then the dataframe is split into \code{Target_data} and \code{Others_data} based on subsetting by the target cells/samples.
#'
#' Finally, a scatter plot is then created using these dataframes.
#' Non-target cells will be displayed in gray if \code{show.others=TRUE},
#' and target cell data is displayed on top, colored and shaped based on the \code{color.var}- and \code{shape.var}-associated data.
#'
#' If \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' If \code{do.hover} is set to \code{TRUE}, the plot is coverted from ggplot to plotly & cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed
#' upon hovering the cursor over the plot.
#'
#' Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item \code{size} and \code{opacity} can be used to adjust the size and transparency of the data points.
#' \item Color can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes can be adjusted with \code{shape.panel}.
#' \item Color and shape labels can be changed using \code{rename.color.groups} and \code{rename.shape.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab-completion lookup.
#' }
#'
#' @seealso
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{x.var}, \code{y.var}, \code{color.var}, \code{shape.var}, and \code{hover.data} options are.
#'
#' \code{\link{dittoDimPlot}} for making very similar data representations, but where dimensionality reduction (PCA, t-SNE, UMAP, etc.) dimensions are the scatterplot axes.
#'
#' @author Daniel Bunis
#' @export
#' @examples
#'
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' # Mock up some percent.mito metadata
#' pbmc$percent.mito <- sample(c(runif(75,0,0.05),runif(5,0.05,0.2)))
#'
#' dittoScatterPlot(
#'     pbmc, x.var = "nCount_RNA", y.var = "nFeature_RNA")
#'
#' # Shapes or colors can be overlaid representing discrete metadata
#' #   or (only colors) continuous metadata / expression data by providing
#' #   metadata or gene names to 'color.var' and 'shape.var'
#' dittoScatterPlot(
#'     pbmc, x.var = "nCount_RNA", y.var = "nFeature_RNA",
#'     color.var = "RNA_snn_res.1",
#'     shape.var = "RNA_snn_res.0.8")
#' dittoScatterPlot(
#'     pbmc, x.var = "nCount_RNA", y.var = "nFeature_RNA",
#'     color.var = "percent.mito",
#'     shape.var = "RNA_snn_res.0.8")
#' dittoScatterPlot(
#'     pbmc, x.var = "nCount_RNA", y.var = "nFeature_RNA",
#'     color.var = "CD14",
#'     shape.var = "RNA_snn_res.0.8")
#'
#' # Note: scatterplots like this can be very useful for dataset QC, especially
#' #   with percentage of reads coming from genes as the color overlay.
dittoScatterPlot <- function(
    object, x.var, y.var, color.var = NULL, shape.var = NULL,
    cells.use = NULL, show.others = FALSE,
    size = 1, opacity = 1,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    assay.x = .default_assay(object), slot.x = .default_slot(object),
    adjustment.x = NULL,
    assay.y = .default_assay(object), slot.y = .default_slot(object),
    adjustment.y = NULL,
    assay.color = .default_assay(object), slot.color = .default_slot(object),
    adjustment.color = NULL,
    do.hover = FALSE, hover.data = NULL, hover.assay = .default_assay(object),
    hover.slot = .default_slot(object), hover.adjustment = NULL,
    shape.panel=c(16,15,17,23,25,8),
    rename.color.groups = NULL, rename.shape.groups = NULL,
    min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
    xlab = x.var, ylab = y.var,
    main = "make", sub = NULL, theme = theme_bw(),
    legend.show = TRUE,
    legend.color.title = color.var, legend.color.size = 5,
    legend.color.breaks = waiver(), legend.color.breaks.labels = waiver(),
    legend.shape.title = shape.var, legend.shape.size = 5,
    data.out = FALSE) {

    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    # Make dataframe
    dat <- data.frame(
        X = .var_OR_get_meta_or_gene(x.var, object, assay.x, slot.x, adjustment.x),
        Y = .var_OR_get_meta_or_gene(y.var, object, assay.y, slot.y, adjustment.y),
        row.names = all.cells)
    aes.args <- list(x = "X", y = "Y")
    do.color <- FALSE
    do.shape <- FALSE
    if (!is.null(color.var)) {
        dat$color <- .var_OR_get_meta_or_gene(
            color.var, object, assay.color, slot.color, adjustment.color)
        dat$color <- .rename_and_or_reorder(
            dat$color, relabels = rename.color.groups)
        aes.args$color = "color"
        do.color <- TRUE
    }
    if (!is.null(shape.var)) {
        dat$shape <- .var_OR_get_meta_or_gene(
            shape.var, object)
        dat$shape <- .rename_and_or_reorder(
            dat$shape, relabels = rename.shape.groups)
        aes.args$shape = "shape"
        do.shape <- TRUE
    }
    if (do.hover) {
        hover.string <- .make_hover_strings_from_vars(
            hover.data, object, hover.assay, hover.slot, hover.adjustment)
        aes.args$text = "hover.string"
    }
    Target_data <- dat[cells.use,]
    Others_data <- dat[!(all.cells %in% cells.use),]

    # Set title
    main <- .leave_default_or_null(main,
        paste0(c(color.var, shape.var), collapse = " and "))

    ### Start building the plot ###
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
    # Add data
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
    if (do.hover) {
        p <- p + suppressWarnings(do.call(geom_point, geom.args))
    } else {
        p <- p + do.call(geom_point, geom.args)
    }
    if (!legend.show) {
        p <- .remove_legend(p)
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
