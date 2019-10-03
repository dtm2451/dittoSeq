################################# dittoPlot ####################################
#'
#' Plots continuous data for cutomizable cells'/samples' groupings on a y-axis
#' @import ggplot2
#'
#' @param var Single string representing the name of a metadata or gene, OR a numeric vector with length equal to the total number of cells/samples in the dataset.
#' This is the data that will be displayed.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups. REQUIRED.
#' @param color.by String representing the name of a metadata to use for setting color.
#' Affects boxplot, vlnplot, and ridgeplot fills.  Default's to \code{group.by} so this input can be skipped if both are the same.
#' @param shape.var Single string representing the name of a metadata to use for setting the shapes of the jitter points.  When not provided, all cells/samples will be represented with dots.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot". See details section for more info.
#' @param data.type String, for when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param do.hover Logical. Default = \code{FALSE}.
#' If set to \code{TRUE} (and if there is a jitter - the data it will work with): object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them,
#' and 'hover.data' argument will be used to determine what data to use.
#'
#' Note: Currently, incompatible with RidgePlots as plotly does not support the geom.
#' @param hover.data String vector, a list of variable names, c("meta1","gene1","meta2","gene2") which determines what data to show upon hover when do.hover is set to \code{TRUE}.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.var}, this sets the panel of shapes which will be used.
#' When nothing is supplied to \code{shape.var}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#' @param main String, sets the plot title. Default = "make" and if left as make, a title will be automatically generated.  To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments. Default = \code{theme_classic()}. See \code{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options.
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param ylab String, sets the continuous-axis label (=y-axis for box and violin plots, x-axis for ridgeplots).
#' Defaults to "\code{var}" or "\code{var} expression" if var is a gene.
#' @param y.breaks Numeric vector, a set of breaks that should be used as major gridlines. c(break1,break2,break3,etc.).
#' @param min,max Scalars which control the zoom of the plot.
#' These inputs set the minimum / maximum values of the data to show.
#' Default = set based on the limits of the data in var.
#' @param x.labels String vector, c("label1","label2","label3",...) which overrides the names of the samples/groups.  NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of x-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of x.reorder should be these indices, but in the order that you would like them rearranged to be.
#' @param x.labels.rotate Logical which sets whether the labels should be rotated.
#' Default: \code{TRUE} for violin and box plots, but \code{FALSE} for ridgeplots.
#' @param add.line numeric value(s) where a dashed horizontal line should go
#' @param line.linetype String which sets the type of line.  Any ggplot linetype should work.  Defaults to "dashed"
#' @param line.color String that sets the color(s) of the horizontal line(s)
#' @param jitter.size Scalar which sets the size of the jitter shapes.
#' @param jitter.width Scalar that sets the width/spread of the jitter in the x direction. Ignored in ridgeplots.
#' @param jitter.color String which sets the color of the jitter shapes
#' @param jitter.shape.legend.size Scalar which changes the size of the shape key in the legend.
#' If set to \code{NA}, \code{jitter.size} is used.
#' @param jitter.shape.legend.show Logical which sets whether the shapes legend will be shown is \code{shape.var} when jitter shape is determined by a variable.
#' @param boxplot.width Scalar which sets the width/spread of the boxplot in the x direction
#' @param boxplot.color String which sets the color of the lines of the boxplot
#' @param boxplot.show.outliers Logical, whether outliers should by including in the boxplot.
#' Default is \code{FALSE} when there is a jitter plotted, \code{TRUE} if there is no jitter.
#' @param boxplot.fill Logical, whether the boxplot should be filled in or not.
#' Note: when boxplot fill is turned off, outliers do not render.
#' @param vlnplot.lineweight Scalar which sets the thickness of the line that outlines the violin plots.
#' @param vlnplot.width Scalar which sets the width/spread of the jitter in the x direction
#' @param vlnplot.scaling String which sets how the widths of the of violin plots are set in relation to eachother.
#' Options are "area", "count", and "width". If the deafult is not right for your data, I recommend trying "width".
#' For a detailed explanation of each, see \code{\link{geom_violin}}.
#' @param ridgeplot.lineweight Scalar which sets the thickness of the ridgeplot outline.
#' @param ridgeplot.scale Scalar which sets the distance/overlap between ridgeplots.
#' A value of 1 means the tallest density curve just touches the baseline of the next higher one.
#' Higher numbers lead to greater overlap.  Default = 1.25
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.title String or \code{NULL}, sets the title for the main legend which includes colors and data representations.
#' This input is set to \code{NULL} by default.
#' @param data.out Logical which sets whether just the plot should be output, or a list containing the plot (\code{p}) and data (\code{data}).  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @param ... arguments passed to dittoPlot by dittoRidgePlot and dittoBoxPlot.  Options are all the ones above.
#' @return a ggplot or plotly where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#' Alternatively, will return the data that would go into such a plot as well with \code{data.out=TRUE}
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided directly, it just uses that) plus metadata associated with the \code{group.by} and \code{color.by} variables.
#' The \code{data.type} input can be used to change what slot of expression data is used when displaying gene expression.
#' If a set of cells to use is indicated with the \code{cells.use} input, the dataframe is then subset to include only those cells.
#' Then, a plot where data is grouped by the \code{group.by} metadata and colored by the \code{color.by} metadata is generated.
#'
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Inclusion of \code{"ridgeplot"} overrides boxplot and violin plot and changes the plot to be horizontal.
#'
#' \code{dittoRidgePlot} and \code{dittoBoxPlot} are included as wrappers of the basic \code{dittoPlot} function.
#' They change the default for the \code{plots} input to be \code{"ridgeplot"} and \code{c("boxplot","jitter")}, respectively.
#'
#' If \code{data.out=TRUE}, a list containing the plot (\code{p}) and a the dataframe containing the underlying data (\code{data}) are returned.
#'
#' If \code{do.hover=TRUE}, cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed
#' upon hovering the cursor over a jitter point.
#'
#' Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item Each data representation has options which are controlled by variables that start with their associated string.
#' For example, all jitter adjustments, like \code{jitter.size}, start with "\code{jitter.}".
#' \item Color can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes used in conjunction with \code{shape.var} can be adjusted with \code{shape.panel}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The legend can be hidden by setting \code{legend.show = TRUE}.
#' \item y-axis zoom and tick marks can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item x-axis labels and groupings can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned off with \code{x.labels.rotate = FALSE}.
#' \item Single or multiple value(s) can be provided to \code{add.line} to add horizontal (or vertical, for ridgeplots) lines at the values provided.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' }
#' @seealso
#' \code{\link{multi_dittoPlot}} for easy creation of multiple dittoPlots that each focus on a different \code{var}
#'
#' \code{\link{dittoPlotVarsAcrossGroups}} to create dittoPlots that show summarized the expression (or value for metadatas), accross groups, of multiple \code{vars} in a single plot.
#'
#' @examples
#' library(Seurat)
#' pbmc <- pbmc_small
#' dittoPlot("CD14", object = "pbmc", group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' dittoPlot("CD14", group.by = "RNA_snn_res.1")
#'
#' # We can adjust the types of plots displayed with the plots input:
#' dittoPlot("CD14", group.by = "RNA_snn_res.1",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     boxplot.fill = FALSE)
#'
#' # Quickly make a Ridgeplot
#' dittoRidgePlot("CD14", group.by = "RNA_snn_res.1")
#'
#' # Quickly make a Boxplot
#' dittoBoxPlot("CD14", group.by = "RNA_snn_res.1")
#'
#' # Any of these can be combined with 'hovering' to retrieve specific info
#' #   about certain data points.  Just add 'do.hover = TRUE' and pick what
#' #   extra data to display by provide set of gene or metadata names to
#' #   'hover.data'.
#' #     Note: ggplotly plots ignores certain dittoSeq plot tweaks.
#' dittoBoxPlot("CD14", group.by = "RNA_snn_res.1",
#'     do.hover = TRUE, hover.data = c("MS4A1","RNA_snn_res.0.8","ident"))
#'
#' @author Daniel Bunis
#' @export

dittoPlot <- function(
    var, object = DEFAULT, group.by, color.by = group.by,
    shape.var = NULL,
    cells.use = NULL, plots = c("jitter","vlnplot"), data.type = "normalized",
    do.hover = FALSE, hover.data = var,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    shape.panel = c(16,15,17,23,25,8),
    theme = theme_classic(), main = "make", sub = NULL,
    ylab = "make", y.breaks = NULL, min = NULL, max = NULL,
    xlab = group.by, x.labels = NULL, x.labels.rotate = NA, x.reorder = NULL,
    jitter.size=1, jitter.width=0.2, jitter.color = "black",
    jitter.shape.legend.size = NA,
    jitter.shape.legend.show = TRUE,
    boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = NA,
    boxplot.fill =TRUE,
    vlnplot.lineweight = 1, vlnplot.width = 1, vlnplot.scaling = "area",
    ridgeplot.lineweight = 1, ridgeplot.scale = 1.25,
    add.line = NULL, line.linetype = "dashed", line.color = "black",
    legend.show = TRUE, legend.title = "make", data.out = FALSE){

    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    #Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    #Establish the full list of cell/sample names
    all.cells <- .all_cells(object)

    #Parse Title Defaults
    exp <- NULL
    if (is.gene(var[1], object)) {
        exp <- " expression"
    }
    ylab <- .leave_default_or_null(ylab,
        default = paste0(var,exp),
        null.if = !(length(var)==1 && is.character(var)))
    main <- .leave_default_or_null(ylab, var,
        null.if = !(length(var)==1 && is.character(var)))
    legend.title <- .leave_default_or_null(legend.title, var,
        null.if = is.null(shape.var))

    #Grab the data
    if (!is.null(shape.var) && is.meta(shape.var, object)) {
        extra.vars = shape.var
    } else {
        extra.vars = NULL
    }
    Target_data <- .dittoPlot_data_gather(
        var, object, group.by, color.by, extra.vars, cells.use, data.type,
        do.hover, hover.data)$Target_data
    Target_data$grouping <-
        .rename_and_or_reorder(Target_data$grouping, x.reorder, x.labels)

    # Make the plot
    p <- ggplot(Target_data, aes_string(fill="color")) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)
    if(!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, Target_data, plots, xlab, ylab, shape.var, jitter.size,
            jitter.width, jitter.color, shape.panel, jitter.shape.legend.size,
            jitter.shape.legend.show, boxplot.width, boxplot.color,
            boxplot.show.outliers, boxplot.fill, vlnplot.lineweight,
            vlnplot.width, vlnplot.scaling, add.line, line.linetype,
            line.color, x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, Target_data, plots, xlab, ylab, jitter.size, jitter.color,
            jitter.shape.legend.size, jitter.shape.legend.show,
            ridgeplot.lineweight, ridgeplot.scale, add.line, line.linetype,
            line.color, x.labels.rotate, do.hover, color.panel, colors,
            y.breaks, min, max)
    }
    #Remove legend, if warrented
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    #DONE. Return the plot or data
    if (data.out) {
        return(list(p = p, data = Target_data))
    } else {
        if (do.hover & ("jitter" %in% plots)) {
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

#### multi_dittoPlot : a function for quickly making multiple DBPlots arranged in a grid.
#' Generates multiple dittoPlots arranged into a grid.
#'
#' @param vars c("var1","var2","var3",...). REQUIRED. A list of vars from which to generate the separate plots
#' @param object the Seurat, SingleCellExperiment, or RNAseq object to draw from, or the "quoted" name of such an object. REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by "metadata" to use for separating values. REQUIRED.
#' @param color.by "metadata" to use for coloring. Affects boxplot, vlnplot, or ridgeplot fills. Defaults to \code{group.by} if not provided.
#' @param show.legend TRUE/FALSE. Whether or not you would like a legend to be plotted.  Default = FALSE
#' @param ncol Integer which sets how many plots will be arranged per row.  Default = 3.
#' @param nrow  Integer which sets how many rows to arrange the plots into.  Default = NULL(/blank) --> becomes however many rows are needed to show all the data.
#' @param add.title Logical which sets whether a title should be added to each individual plot
#' @param ylab Logical or string (with "var" being a special case).
#' If logical, sets whether a y-axis title should be added.
#' When set to \code{TRUE}, default dittoPlot behavior will be observed: y-label for any gene vars will be "'var' expression"
#' When set to \code{"var"}, then the \code{vars} names alone will be used.
#' When set as any other string, that string will be used as the y-axis label for every plot.
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{NULL}, which removes it from plotting.
#' @param OUT.List Logical which sets whether the output should be a list of objects instead of the plots arranged into a single plot grid.
#' Outputting as list allows manual input into gridArrange for moving plots around / adjusting sizes.
#' In the list, all plots will be named by the element of \code{vars} the represent.
#' @param ... other paramters passed along to dittoPlot.
#' @return Given multiple 'var' parameters, this function will output a DBPlot for each one, arranged into a grid.  All parameters that can be adjusted in DBPlot can be adjusted here.
#' @seealso
#' \code{\link{dittoPlot}} for the single plot version of this function
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' genes <- c("CD8A","CD3E","FCER1A","CD14")
#' multi_dittoPlot(genes, object = "pbmc",
#'     group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object
#'       # input can be skipped completely.
#' DEFAULT <- "pbmc"
#' multi_dittoPlot(genes,
#'     group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#'
#' #To make it output a grid that is 2x2, to add y-axis labels
#' # instead of titles, and to show legends...
#' multi_dittoPlot(genes,
#'     group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
#'     nrow = 2, ncol = 2,
#'     add.title = FALSE, ylab = TRUE,  #Add y axis labels instead of titles
#'     show.legend = TRUE)              #Show legends
#' # To eliminate the "expression", change ylab = TRUE to ylab = "var"
#' multi_dittoPlot(genes,
#'             group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
#'             nrow = 2, ncol = 2,              #Make it 2x2
#'             add.title = FALSE, ylab = "var", #Add y axis labels without "expression"
#'             show.legend = TRUE)              #Show legends
#'
#' @author Daniel Bunis
#' @export

multi_dittoPlot <- function(
    vars, object = DEFAULT, group.by, color.by = group.by, show.legend = FALSE,
    ncol = 3, nrow = NULL, add.title=TRUE, ylab = FALSE, xlab = NULL, OUT.List = FALSE,
    ...) {

    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }

    ylab.input <- ylab

    plots <- lapply(vars, function(X) {
        dittoPlot(X, object, group.by, color.by, xlab = xlab,
            ylab =
                if (is.character(ylab.input)) {
                    ifelse(
                        ylab.input == "var",
                        X,
                        ylab.input)
                } else {
                    if (ylab.input) {
                        "make" }
                    else {
                        NULL
                    }
                },
            main =
                if (add.title) {
                    "make"
                } else {
                    NULL
                },
            ...) +
            theme(legend.position = ifelse(show.legend, "right", "none"))
    })

    #Output
    if (OUT.List){
        names(plots) <- vars
        return(plots)
    } else {
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings horizontally in a desnsity representation
#' @export
dittoRidgePlot <- function(..., plots = c("ridgeplot")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings in boxplot form
#' @export
dittoBoxPlot <- function(..., plots = c("boxplot","jitter")){ dittoPlot(..., plots = plots) }

.dittoPlot_add_data_y_direction <- function(
    p, Target_data, plots, xlab, ylab, shape.var,
    jitter.size, jitter.width, jitter.color,shape.panel,
    jitter.shape.legend.size, jitter.shape.legend.show,
    boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
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
                scale = vlnplot.scaling)
        }
        if (plots[i] == "boxplot") {
            boxplot.args <- list(
                width=boxplot.width,
                color = boxplot.color,
                alpha = ifelse(boxplot.fill, 1, 0))
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
            #If shape.var metadata given, use it. Else, shapes[1] which = dots (16) by default
            if (!is.null(shape.var) && is.meta(shape.var, object)) {
                #Make jitter with shapes
                jitter.args$mapping <- if (do.hover) {
                    aes_string(shape = shape.var, text = "hover.string")
                } else {
                    aes_string(shape = shape.var)
                }
                p <- p + do.call(geom_jitter, jitter.args) +
                    scale_shape_manual(
                        values = shape.panel[seq_along(meta.levels(
                            shape.var, object, rownames(Target_data)))])
                if (!is.na(jitter.shape.legend.size)){
                    p <- p + guides(shape = guide_legend(
                        override.aes = list(size=jitter.shape.legend.size)))
                }
                if (jitter.shape.legend.show==FALSE){
                    p <- p + guides(shape = "none")
                }
            } else {
                if (do.hover) {
                    jitter.args$mapping <- aes_string(text = "hover.string")
                }
                jitter.args$shape <- shape.panel[1]
                p <- p + do.call(geom_jitter, jitter.args)
            }
        }
    }

    # Add labels and, if requested, lines
    p <- p + aes_string(x = "grouping", y = "var.data") +
        xlab(xlab) + ylab(ylab)
    if (is.na(x.labels.rotate) || x.labels.rotate) {
        p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1, size=12))
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
                            add.line=NULL, line.linetype = "dashed", line.color = "black",
                            x.labels.rotate = FALSE, do.hover, color.panel, colors, y.breaks, min, max) {
    #This function takes in a partial dittoPlot ggplot object without any data overlay, and parses adding the main data visualizations.
    # It adds plots based on what is requested in plots, *ordered by their order*

    # Now that we know the plot's direction, set y-axis limits
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
    p <- p + aes_string(x = "var.data", y = "grouping") + xlab(ylab) + ylab(xlab) +
        scale_y_discrete(expand = expand_scale(mult=c(0, 0.68)))
    if (!is.na(x.labels.rotate) && x.labels.rotate) {
        p <- p + theme(axis.text.y= element_text(angle=45, hjust = 1, vjust = 1, size=12))
    }
    if (!is.null(add.line)) {
        p <- p + geom_vline(xintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

.dittoPlot_data_gather <- function(
    main.var, object = DEFAULT, group.by = "Sample", color.by = group.by,
    extra.vars = NULL, cells.use = NULL, data.type = "normalized",
    do.hover = FALSE, hover.data = c(main.var, extra.vars)) {

    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    # Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    # Establish the full list of cell/sample names
    all.cells <- .all_cells(object)
    ### Make dataframe for storing the plotting data:
    full_data <- data.frame(
        var.data = .var_OR_get_meta_or_gene(main.var, object, data.type),
        grouping = meta(group.by, object),
        color = meta(color.by, object),
        row.names = all.cells)
    names <- names(full_data)
    # Add Extra data
    if(length(extra.vars)>0){
        for (i in seq_along(extra.vars)){
            full_data <- cbind(full_data, .var_OR_get_meta_or_gene(extra.vars, object, data.type))
        }
        names <- c(names, extra.vars)
    }
    # Add hover strings
    if (do.hover) {
        full_data$hover.string <- .make_hover_strings_from_vars(hover.data, object, data.type)
        names <- c(names, "hover.string")
    }
    # Add column names
    colnames(full_data) <- names

    return(list(Target_data = full_data[all.cells %in% cells.use,],
                Others_data = full_data[!(all.cells %in% cells.use),]))
}

