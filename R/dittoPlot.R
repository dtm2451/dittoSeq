################################# dittoPlot ####################################
#'
#' Plots continuous data for cutomizable cells'/samples' groupings on a y-axis
#' @import ggplot2
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param var Single string representing the name of a metadata or gene, OR a vector with length equal to the total number of cells/samples in the dataset.
#' This is the data that will be displayed.
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param color.by String representing the name of a metadata to use for setting fills.
#' Great for highlighting subgroups when wanted, but it defaults to \code{group.by} so this input can be skipped otherwise.
#' Affects boxplot, vlnplot, and ridgeplot fills.
#' @param shape.by Single string representing the name of a metadata to use for setting the shapes of the jitter points.  When not provided, all cells/samples will be represented with dots.
#' @param split.by 1 or 2 strings naming discrete metadata to use for splitting the cells/samples into multiple plots with ggplot faceting.
#'
#' When 2 metadatas are named, c(row,col), the first is used as rows and the second is used for columns of the resulting grid.
#'
#' When 1 metadata is named, shape control can be achieved with \code{split.nrow} and \code{split.ncol}
#'
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}.
#' @param extra.vars String vector providing names of any extra metadata to be stashed in the dataframe supplied to \code{ggplot(data)}.
#'
#' Useful for making custom spliting/faceting or other additional alterations \emph{after} dittoSeq plot generation.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot".
#' Order matters: c("vlnplot", "boxplot", "jitter") will put a violin plot in the back, boxplot in the middle, and then individual dots in the front.
#' See details section for more info.
#' @param assay,slot single strings or integer that set which data to use when plotting gene expression / feature data. See \code{\link{gene}} for more information.
#' @param adjustment When plotting gene expression / feature counts, should that data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param swap.rownames String. For SummarizeedExperiment or SingleCellExperiment objects, the column name of rowData(object) to be used to identify features instead of rownames(object).
#' @param do.hover Logical. Default = \code{FALSE}.
#' If set to \code{TRUE} (and if there is a "jitter" in \code{plots}): object will be converted to a ggplotly object so that data about individual cells will be displayed when you hover your cursor over the jitter points,
#'
#' Note: Currently, hovering is incompatible with RidgePlots as plotly does not support the geom_density_ridges2 geom.
#' @param hover.data String vector, a list of variable names, c("meta1","gene1","meta2",...) which determines what data to show upon hover when do.hover is set to \code{TRUE}.
#' @param color.panel String vector which sets the colors to draw from for plot fills.
#' Default = \code{dittoColors()}.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use.
#' (Provides an alternative to directly modifying \code{color.panel}.)
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.by}, this sets the panel of shapes which will be used.
#' When nothing is supplied to \code{shape.by}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#' @param main String, sets the plot title. Default = "make" and if left as make, a title will be automatically generated.  To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_classic()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Set to \code{NULL} to remove.
#' @param ylab String, sets the continuous-axis label (=y-axis for box and violin plots, x-axis for ridgeplots).
#' Defaults to "\code{var}" or "\code{var} expression" if \code{var} is a gene.
#' @param y.breaks Numeric vector, a set of breaks that should be used as major gridlines. c(break1,break2,break3,etc.).
#' @param min,max Scalars which control the zoom of the plot.
#' These inputs set the minimum / maximum values of the data to show.
#' Default = set based on the limits of the data in var.
#' @param x.labels String vector, c("label1","label2","label3",...) which overrides the names of the samples/groups.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of x-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of x.reorder should be these indices, but in the order that you would like them rearranged to be.
#' 
#' Recommendation for advanced users: If you find yourself coming back to this input too many times, an alternative solution that can be easier long-term
#' is to make the target data into a factor, and to put its levels in the desired order: \code{factor(data, levels = c("level1", "level2", ...))}.
#' \code{\link{metaLevels}} can be used to quickly get the identities that need to be part of this 'levels' input.
#' @param x.labels.rotate Logical which sets whether the labels should be rotated.
#' Default: \code{TRUE} for violin and box plots, but \code{FALSE} for ridgeplots.
#' @param add.line numeric value(s) where one or multiple line should be added
#' @param line.linetype String which sets the type of line for \code{add.line}.
#' Defaults to "dashed", but any ggplot linetype will work.
#' @param line.color String that sets the color(s) of the \code{add.line} line(s)
#' @param jitter.size Scalar which sets the size of the jitter shapes.
#' @param jitter.width Scalar that sets the width/spread of the jitter in the x direction. Ignored in ridgeplots.
#' @param jitter.color String which sets the color of the jitter shapes
#' @param jitter.shape.legend.size Scalar which changes the size of the shape key in the legend.
#' If set to \code{NA}, \code{jitter.size} is used.
#' @param jitter.shape.legend.show Logical which sets whether the shapes legend will be shown when its shape is determined by \code{shape.by}.
#' @param do.raster Logical. When set to \code{TRUE}, rasterizes the jitter plot layer, changing it from individually encoded points to a flattened set of pixels.
#' This can be useful for editing in external programs (e.g. Illustrator) when there are many thousands of data points.
#' @param raster.dpi Number indicating dots/pixels per inch (dpi) to use for rasterization. Default = 300.
#' @param boxplot.width Scalar which sets the width/spread of the boxplot in the x direction
#' @param boxplot.color String which sets the color of the lines of the boxplot
#' @param boxplot.show.outliers Logical, whether outliers should by including in the boxplot.
#' Default is \code{FALSE} when there is a jitter plotted, \code{TRUE} if there is no jitter.
#' @param boxplot.fill Logical, whether the boxplot should be filled in or not.
#' Known bug: when boxplot fill is turned off, outliers do not render.
#' @param boxplot.position.dodge Scalar which adjusts the distance between boxplots when multiple are drawn per grouping (a.k.a. when \code{group.by} and \code{color.by} are not equal).
#' @param vlnplot.lineweight Scalar which sets the thickness of the line that outlines the violin plots.
#' @param vlnplot.width Scalar which sets the width/spread of the jitter in the x direction
#' @param vlnplot.scaling String which sets how the widths of the of violin plots are set in relation to eachother.
#' Options are "area", "count", and "width". If the deafult is not right for your data, I recommend trying "width".
#' For a detailed explanation of each, see \code{\link{geom_violin}}.
#' @param ridgeplot.lineweight Scalar which sets the thickness of the ridgeplot outline.
#' @param ridgeplot.scale Scalar which sets the distance/overlap between ridgeplots.
#' A value of 1 means the tallest density curve just touches the baseline of the next higher one.
#' Higher numbers lead to greater overlap.  Default = 1.25
#' @param ridgeplot.ymax.expansion Scalar which adjusts the minimal space between the topmost grouping and the top of the plot in order to ensure the curve is not cut off by the plotting grid.
#' The larger the value, the greater the space requested.
#' When left as NA, dittoSeq will attempt to determine an ideal value itself based on the number of groups & linear interpolation between these goal posts: #groups of 3 or fewer: 0.6; #groups=12: 0.1; #groups or 34 or greater: 0.05.
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.title String or \code{NULL}, sets the title for the main legend which includes colors and data representations.
#' This input is set to \code{NULL} by default.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot (\code{p}) and data (\code{data}).
#'
#' Note: plotly conversion is turned off in the \code{data.out = TRUE} setting, but hover.data is still calculated.
#' @param ... arguments passed to dittoPlot by dittoRidgePlot, dittoRidgeJitter, and dittoBoxPlot wrappers.
#' Options are all the ones above.
#' @return a ggplot or plotly where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively when \code{data.out=TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the plot where additional data will be displayed when the cursor is hovered over jitter points.
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided, that data).
#' On the discrete axis, data will be grouped by the metadata given to \code{group.by} and colored by the metadata given to \code{color.by}.
#' The \code{assay} and \code{slot} inputs can be used to change what expression data is used when displaying gene expression.
#' If a set of cells to use is indicated with the \code{cells.use} input, the data is subset to include only those cells before plotting.
#'
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Inclusion of \code{"ridgeplot"} overrides boxplot and violin plot and changes the plot to be horizontal.
#'
#' When \code{split.by} is provided the name of a metadata containing discrete data, separate plots will be produced representing each of the distinct groupings of the split.by data.
#'
#' \code{dittoRidgePlot}, \code{dittoRidgeJitter}, and \code{dittoBoxPlot} are included as wrappers of the basic \code{dittoPlot} function
#' that simply change the default for the \code{plots} input to be \code{"ridgeplot"}, \code{c("ridgeplot","jitter")}, or \code{c("boxplot","jitter")},
#' to make such plots even easier to produce.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item Each data representation has options which are controlled by variables that start with their associated string.
#' For example, all jitter adjustments, like \code{jitter.size}, start with "\code{jitter.}".
#' \item Colors can be adjusted with \code{color.panel}.
#' \item Shapes used in conjunction with \code{shape.by} can be adjusted with \code{shape.panel}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The legend can be hidden by setting \code{legend.show = TRUE}.
#' \item y-axis zoom and tick marks can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item x-axis labels and groupings can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned off with \code{x.labels.rotate = FALSE}.
#' \item Line(s) can be added at single or multiple value(s) by providing these values to \code{add.line}.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' \item Single or multiple additional per-cell features can be retrieved and stashed within the underlying data using \code{extra.vars}.
#' This can be very useful for making manual additional alterations \emph{after} dittoSeq plot generation.
#' }
#' @seealso
#' \code{\link{multi_dittoPlot}} for easy creation of multiple dittoPlots each focusing on a different \code{var}.
#'
#' \code{\link{dittoPlotVarsAcrossGroups}} to create dittoPlots that show summarized expression (or values for metadata), accross groups, of multiple \code{vars} in a single plot.
#'
#' \code{\link{dittoRidgePlot}}, \code{\link{dittoRidgeJitter}}, and \code{\link{dittoBoxPlot}} for shortcuts to a few 'plots' input shortcuts
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' # Basic dittoplot, with jitter behind a vlnplot (looks better with more cells)
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint")
#'
#' # Color distinctly from the grouping variable using 'color.by'
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     color.by = "conditions")
#'
#' # Update the 'plots' input to change / reorder the data representations
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"))
#'
#' # Modify the look with intuitive inputs
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     boxplot.color = "white",
#'     main = "CD3E",
#'     legend.show = FALSE)
#'
#' # Data can also be split in other ways with 'shape.by' or 'split.by'
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     shape.by = "clustering",
#'     split.by = "SNP") # single split.by element
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     split.by = c("groups","SNP")) # row and col split.by elements
#'
#' # For faceting, instead of using 'split.by', the target data can alternatively
#' # be given to 'extra.var' to have it added in the underlying dataframe, then
#' # faceting can be added manually for extra flexibility
#' dittoPlot(myRNA, "gene1", "clustering",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     extra.var = "SNP") + facet_wrap("SNP", ncol = 1, strip.position = "left")
#'
#' # Quickly make a Ridgeplot
#' dittoRidgePlot(myRNA, "gene1", group.by = "timepoint")
#'
#' # Quickly make a Boxplot
#' dittoBoxPlot(myRNA, "gene1", group.by = "timepoint")
#'
#' @author Daniel Bunis
#' @export

dittoPlot <- function(
    object,
    var,
    group.by,
    color.by = group.by,
    shape.by = NULL,
    split.by = NULL,
    extra.vars = NULL,
    cells.use = NULL,
    plots = c("jitter","vlnplot"),
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL,
    do.hover = FALSE,
    hover.data = var,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    shape.panel = c(16,15,17,23,25,8),
    theme = theme_classic(),
    main = "make",
    sub = NULL,
    ylab = "make",
    y.breaks = NULL,
    min = NULL,
    max = NULL,
    xlab = group.by,
    x.labels = NULL,
    x.labels.rotate = NA,
    x.reorder = NULL,
    split.nrow = NULL,
    split.ncol = NULL,
    do.raster = FALSE,
    raster.dpi = 300,
    jitter.size = 1,
    jitter.width = 0.2,
    jitter.color = "black",
    jitter.shape.legend.size = NA,
    jitter.shape.legend.show = TRUE,
    boxplot.width = 0.2,
    boxplot.color = "black",
    boxplot.show.outliers = NA,
    boxplot.fill = TRUE,
    boxplot.position.dodge = vlnplot.width,
    vlnplot.lineweight = 1,
    vlnplot.width = 1,
    vlnplot.scaling = "area",
    ridgeplot.lineweight = 1,
    ridgeplot.scale = 1.25,
    ridgeplot.ymax.expansion = NA,
    add.line = NULL,
    line.linetype = "dashed",
    line.color = "black",
    legend.show = TRUE,
    legend.title = "make",
    data.out = FALSE) {

    #Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    #Establish the full list of cell/sample names
    all.cells <- .all_cells(object)

    #Parse Title Defaults
    exp <- NULL
    if (isGene(var[1], object, assay)) {
        exp <- " expression"
    }
    ylab <- .leave_default_or_null(ylab,
        default = paste0(var,exp),
        null.if = !(length(var)==1 && is.character(var)))
    main <- .leave_default_or_null(main, var,
        null.if = !(length(var)==1 && is.character(var)))
    legend.title <- .leave_default_or_null(legend.title, var,
        null.if = is.null(shape.by))

    # Grab the data
    Target_data <- .dittoPlot_data_gather(object, var, group.by, color.by,
        c(shape.by,split.by,extra.vars), cells.use, assay, slot, adjustment,
        swap.rownames, do.hover, hover.data)$Target_data
    Target_data$grouping <-
        .rename_and_or_reorder(Target_data$grouping, x.reorder, x.labels)

    # Make the plot
    p <- ggplot(Target_data, aes_string(fill="color")) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)
    if(!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, Target_data, plots, xlab, ylab, shape.by, jitter.size,
            jitter.width, jitter.color, shape.panel, jitter.shape.legend.size,
            jitter.shape.legend.show, do.raster, raster.dpi,
            boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
            boxplot.position.dodge,
            vlnplot.lineweight, vlnplot.width, vlnplot.scaling,
            add.line, line.linetype, line.color,
            x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, Target_data, plots, xlab, ylab, jitter.size, jitter.color,
            jitter.shape.legend.size, jitter.shape.legend.show,
            ridgeplot.lineweight, ridgeplot.scale, ridgeplot.ymax.expansion,
            add.line, line.linetype,
            line.color, x.labels.rotate, do.hover, color.panel, colors,
            y.breaks, min, max)
    }
    # Extra tweaks
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, object, cells.use)
    }
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    #DONE. Return the plot or data
    if (data.out) {
        return(list(p = p, data = Target_data))
    } else {
        if (do.hover & ("jitter" %in% plots)) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings horizontally in a density representation
#' @export
dittoRidgePlot <- function(..., plots = c("ridgeplot")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot dittoRidgePlot, but with jitter overlaid
#' @export
dittoRidgeJitter <- function(..., plots = c("ridgeplot", "jitter")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings in boxplot form
#' @export
dittoBoxPlot <- function(..., plots = c("boxplot","jitter")){ dittoPlot(..., plots = plots) }

.dittoPlot_add_data_y_direction <- function(
    p, Target_data, plots, xlab, ylab, shape.by,
    jitter.size, jitter.width, jitter.color,shape.panel,
    jitter.shape.legend.size, jitter.shape.legend.show, do.raster, raster.dpi,
    boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
    boxplot.position.dodge,
    vlnplot.lineweight, vlnplot.width, vlnplot.scaling, add.line,
    line.linetype, line.color, x.labels.rotate, do.hover, y.breaks, min, max,
    object) {
    # This function takes in a partial dittoPlot ggplot object without any data
    # overlay, and parses adding the main data visualizations.
    # Adds plots based on what is requested in plots, ordered by their order.

    # Now that we know the plot's direction, set y-axis limits
    if (!is.null(y.breaks)) {
        p <- p + scale_y_continuous(breaks = y.breaks)
    }
    if (is.null(min)) {
        min <- min(Target_data$var.data)
    }
    if (is.null(max)) {
        max <- max(Target_data$var.data)
    }
    p <- p + coord_cartesian(ylim=c(min,max))

    # Add Plots
    for (i in seq_along(plots)) {
        if (plots[i] == "vlnplot") {
            p <- p + geom_violin(
                size = vlnplot.lineweight,
                width = vlnplot.width,
                scale = vlnplot.scaling,
                na.rm = TRUE)
        }

        if (plots[i] == "boxplot") {
            boxplot.args <- list(
                width = boxplot.width,
                color = boxplot.color,
                alpha = ifelse(boxplot.fill, 1, 0),
                position = position_dodge(width = boxplot.position.dodge),
                na.rm = TRUE)
            if (is.na(boxplot.show.outliers)) {
                boxplot.show.outliers <- ifelse("jitter" %in% plots, FALSE, TRUE)
            }
            if (!boxplot.show.outliers) {
                boxplot.args$outlier.shape <- NA
            }
            p <- p + do.call(geom_boxplot, boxplot.args)
        }

        if (plots[i] == "jitter") {
            jitter.args <- list(
                size=jitter.size,
                width=jitter.width,
                height = 0,
                color = jitter.color)
            
            geom_for_jitter <- geom_jitter
            if (do.raster) {
                .error_if_no_ggrastr()
                geom_for_jitter <- ggrastr::geom_jitter_rast
                jitter.args$raster.dpi <- raster.dpi
            }
            
            jitter.aes.args <- list()
            if (do.hover) {
                jitter.aes.args$text <-  "hover.string"
            }
            
            #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
            if (!is.null(shape.by) && isMeta(shape.by, object)) {
                
                # Set shape in aes & set scales/theming.
                jitter.aes.args$shape <- shape.by
                
                p <- p + scale_shape_manual(
                    values = shape.panel[seq_along(metaLevels(shape.by, object, rownames(Target_data)))])
                
                if (!is.na(jitter.shape.legend.size)){
                    p <- p + guides(shape = guide_legend(
                        override.aes = list(size=jitter.shape.legend.size)))
                }
                if (jitter.shape.legend.show==FALSE){
                    p <- p + guides(shape = "none")
                }
                
            } else {
                # Set shape outside of aes
                jitter.args$shape <- shape.panel[1]
            }
            
            jitter.args$mapping <- do.call(aes_string, jitter.aes.args)
            
            if (do.hover) {
                p <- p + suppressWarnings(do.call(geom_for_jitter, jitter.args))
            } else {
                p <- p + do.call(geom_for_jitter, jitter.args)
            }
        }
    }

    # Add labels and, if requested, lines
    p <- p + aes_string(x = "grouping", y = "var.data") +
        xlab(xlab) + ylab(ylab)
    if (is.na(x.labels.rotate) || x.labels.rotate) {
        p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1))
    }
    if (!is.null(add.line)) {
        p <- p + geom_hline(yintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

#' @importFrom ggridges geom_density_ridges2
.dittoPlot_add_data_x_direction <- function(p, Target_data, plots, xlab, ylab,
                            jitter.size=1, jitter.color = "black",
                            jitter.shape.legend.size, jitter.shape.legend.show,
                            ridgeplot.lineweight = 1, ridgeplot.scale = 1.25,
                            ridgeplot.ymax.expansion = NA,
                            add.line=NULL, line.linetype = "dashed", line.color = "black",
                            x.labels.rotate = FALSE, do.hover, color.panel, colors, y.breaks, min, max) {
    #This function takes in a partial dittoPlot ggplot object without any data overlay, and parses adding the main data visualizations.
    # It adds plots based on what is requested in plots, *ordered by their order*

    # Now that we know the plot's direction, set "y"-axis limits
    if (!is.null(y.breaks)) {
        p <- p + scale_x_continuous(breaks = y.breaks)
    }
    if (is.null(min)) {
        min <- min(Target_data$var.data)
    }
    if (is.null(max)) {
        max <- max(Target_data$var.data)
    }
    p <- p + coord_cartesian(xlim=c(min,max))
    
    # For stylistic issues with plotting defaults, also adjust grouping-axis limits
    if (is.na(ridgeplot.ymax.expansion)) {
        num_groups <- length(unique(Target_data$grouping))
        # From 0.6 to 0.1 between 4 to 12 groups and 0.05 by 34 groups.
        set_exp <- stats::approxfun(
            x=c(3,12,34), y=c(0.6, 0.1, 0.05), yleft = 0.6, yright = 0.05)
        ridgeplot.ymax.expansion <- set_exp(num_groups)
    }
    p <- p + 
        scale_y_discrete(expand = expansion(mult=c(0, ridgeplot.ymax.expansion)))

    # Add ridgeplot (and jitter) data
    p <- p + ggridges::geom_density_ridges2(
        size = ridgeplot.lineweight, scale = ridgeplot.scale,
        jittered_points = "jitter" %in% plots, point_size = jitter.size, point_color = jitter.color) +
        scale_color_manual(values=color.panel[colors])
    if (!is.na(jitter.shape.legend.size)) {
        p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
    }
    if (jitter.shape.legend.show==FALSE){
        p <- p + guides(shape = "none")
    }

    # Add labels and, if requested, lines
    p <- p + aes_string(x = "var.data", y = "grouping") + xlab(ylab) + ylab(xlab)
    if (!is.na(x.labels.rotate) && x.labels.rotate) {
        p <- p + theme(axis.text.y= element_text(angle=45, hjust = 1, vjust = 1))
    }
    if (!is.null(add.line)) {
        p <- p + geom_vline(xintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

.dittoPlot_data_gather <- function(
    object, main.var, group.by = "Sample", color.by = group.by,
    extra.vars = NULL, cells.use = NULL,
    assay, slot, adjustment,
    swap.rownames = NULL,
    do.hover = FALSE, hover.data = c(main.var, extra.vars)) {

    # Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    # Establish the full list of cell/sample names
    all.cells <- .all_cells(object)
    ### Make dataframe for storing the plotting data:
    full_data <- data.frame(
        var.data = .var_OR_get_meta_or_gene(
            main.var, object, assay, slot, adjustment, swap.rownames),
        grouping = meta(group.by, object),
        color = meta(color.by, object),
        row.names = all.cells)
    # Add split and extra data
    full_data <- .add_by_cell(full_data, extra.vars, extra.vars, object, assay,
        slot, adjustment, mult = TRUE)
    # Add hover strings
    if (do.hover) {
        full_data$hover.string <- .make_hover_strings_from_vars(
            hover.data, object, assay, slot, adjustment)
    }

    return(list(Target_data = full_data[all.cells %in% cells.use,],
                Others_data = full_data[!(all.cells %in% cells.use),]))
}

