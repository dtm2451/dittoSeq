#' Plots continuous data for customizeable cells'/samples' groupings on a y- (or x-) axis
#' @import ggplot2
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param var Single string representing the name of a metadata or gene, OR a vector with length equal to the total number of cells/samples in the dataset.
#' Alternatively, a string vector naming multiple genes or metadata.
#' This is the primary data that will be displayed.
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param color.by String representing the name of a metadata to use for setting fills.
#' Great for highlighting supersets or subgroups when wanted, but it defaults to \code{group.by} so this input can be skipped otherwise.
#' @param shape.by Single string representing the name of a metadata to use for setting the shapes of the jitter points.  When not provided, all cells/samples will be represented with dots.
#' @param multivar.aes "split", "group", or "color", the plot feature to utilize for displaying 'var' value when \code{var} is given multiple genes or metadata.
#' When set to "split", inputs \code{split.nrow}, \code{split.ncol}, and \code{split.adjust} can be used to 
#' @param multivar.split.dir "row" or "col", sets the direction of faceting used for 'var' values when \code{var} is given multiple genes or metadata, when \code{multivar.aes = "split"}, and when \code{split.by} is used to provide additional data to facet by.
#' @param split.by 1 or 2 strings naming discrete metadata to use for splitting the cells/samples into multiple plots with ggplot faceting.
#'
#' When 2 metadatas are named, c(row,col), the first is used as rows and the second is used for columns of the resulting grid.
#'
#' When 1 metadata is named, shape control can be achieved with \code{split.nrow} and \code{split.ncol}
#'
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}, or when multiple genes/metadata are given to \code{var} and \code{multivar.aes = "split"}.
#' @param split.adjust A named list which allows extra parameters to be pushed through to the faceting function call.
#' List elements should be valid inputs to the faceting functions, e.g. `list(scales = "free")`.
#' 
#' For options, when giving 1 metadata to \code{split.by} or if faceting is by a set of \code{var}s, see \code{\link[ggplot2]{facet_wrap}},
#' OR when giving 2 metadatas to \code{split.by}, see \code{\link[ggplot2]{facet_grid}}.
#' @param extra.vars String vector providing names of any extra metadata to be stashed in the dataframe supplied to \code{ggplot(data)}.
#'
#' Useful for making custom spliting/faceting or other additional alterations \emph{after} dittoSeq plot generation.
#' @param cells.use String vector of cells'(/samples' for bulk data) names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot".
#' 
#' Order matters: c("vlnplot", "boxplot", "jitter") will put a violin plot in the back, boxplot in the middle, and then individual dots in the front.
#' 
#' See details section for more info.
#' @param assay,slot single strings or integers (SCEs and SEs) or an optionally named vector of such values that set which expression data to use.
#' See \code{\link{GeneTargeting}} for specifics and examples -- Seurat and SingleCellExperiment objects deal with these differently, and functionality additions in dittoSeq have led to some minimal divergence from the native methodologies.
#' @param adjustment When plotting gene expression / feature counts, should that data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param do.hover Logical. Default = \code{FALSE}.
#' If set to \code{TRUE}: object will be converted to a ggplotly object so that data about individual cells will be displayed when you hover your cursor over the jitter points (assuming that there is a "jitter" in \code{plots}),
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
#' These inputs set the minimum / maximum values of the data to display.
#' Default = NA, which allows ggplot to set these limits based on the range of all data being shown.
#' @param x.labels String vector, c("label1","label2","label3",...) which overrides the names of groupings.
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
#' @param add.line numeric value(s) where one or multiple line(s) should be added
#' @param line.linetype String which sets the type of line for \code{add.line}.
#' Defaults to "dashed", but any ggplot linetype will work.
#' @param line.color String that sets the color(s) of the \code{add.line} line(s)
#' @param jitter.size Scalar which sets the size of the jitter shapes.
#' @param jitter.width Scalar that sets the width/spread of the jitter in the x direction. Ignored in ridgeplots.
#' 
#' Note for when \code{color.by} is used to split x-axis groupings into additional bins: ggplot does not shrink jitter widths accordingly, so be sure to do so yourself!
#' Ideally, needs to be 0.5/num_subgroups.
#' @param jitter.color String which sets the color of the jitter shapes
#' @param jitter.shape.legend.size Scalar which changes the size of the shape key in the legend.
#' If set to \code{NA}, \code{jitter.size} is used.
#' @param jitter.shape.legend.show Logical which sets whether the shapes legend will be shown when its shape is determined by \code{shape.by}.
#' @param jitter.position.dodge Scalar which adjusts the relative distance between jitter widths when multiple subgroups exist per \code{group.by} grouping (a.k.a. when \code{group.by} and \code{color.by} are not equal).
#' Similar to \code{boxplot.position.dodge} input & defaults to the value of that input so that BOTH will actually be adjusted when only, say, \code{boxplot.position.dodge = 0.3} is given.
#' @param do.raster Logical. When set to \code{TRUE}, rasterizes the jitter plot layer, changing it from individually encoded points to a flattened set of pixels.
#' @param do.raster Logical. When set to \code{TRUE}, rasterizes the jitter plot layer, changing it from individually encoded points to a flattened set of pixels.
#' This can be useful for editing in external programs (e.g. Illustrator) when there are many thousands of data points.
#' @param raster.dpi Number indicating dots/pixels per inch (dpi) to use for rasterization. Default = 300.
#' @param boxplot.width Scalar which sets the width/spread of the boxplot in the x direction
#' @param boxplot.color String which sets the color of the lines of the boxplot
#' @param boxplot.show.outliers Logical, whether outliers should by including in the boxplot.
#' Default is \code{FALSE} when there is a jitter plotted, \code{TRUE} if there is no jitter.
#' @param boxplot.outlier.size Scalar which adjusts the size of points used to mark outliers
#' @param boxplot.fill Logical, whether the boxplot should be filled in or not.
#' Known bug: when boxplot fill is turned off, outliers do not render.
#' @param boxplot.position.dodge Scalar which adjusts the relative distance between boxplots when multiple are drawn per grouping (a.k.a. when \code{group.by} and \code{color.by} are not equal).
#' By default, this input actually controls the value of \code{jitter.position.dodge} unless the \code{jitter} version is provided separately.
#' @param boxplot.lineweight Scalar which adjusts the thickness of boxplot lines.
#' @param vlnplot.lineweight Scalar which sets the thickness of the line that outlines the violin plots.
#' @param vlnplot.width Scalar which sets the width/spread of violin plots in the x direction
#' @param vlnplot.scaling String which sets how the widths of the of violin plots are set in relation to each other.
#' Options are "area", "count", and "width". If the default is not right for your data, I recommend trying "width".
#' For an explanation of each, see \code{\link{geom_violin}}.
#' @param vlnplot.quantiles Single number or numeric vector of values in [0,1] naming quantiles at which to draw a horizontal line within each violin plot. Example: \code{c(0.1, 0.5, 0.9)}
#' @param ridgeplot.lineweight Scalar which sets the thickness of the ridgeplot outline.
#' @param ridgeplot.scale Scalar which sets the distance/overlap between ridgeplots.
#' A value of 1 means the tallest density curve just touches the baseline of the next higher one.
#' Higher numbers lead to greater overlap.  Default = 1.25
#' @param ridgeplot.ymax.expansion Scalar which adjusts the minimal space between the top-most grouping and the top of the plot in order to ensure that the curve is not cut off by the plotting grid.
#' The larger the value, the greater the space requested.
#' When left as NA, dittoSeq will attempt to determine an ideal value itself based on the number of groups & linear interpolation between these goal posts: 0.6 when g<=3, 0.1 when g==12, and 0.05 when g>=34, where g is the number of groups.
#' @param ridgeplot.shape Either "smooth" or "hist", sets whether ridges will be smoothed (the typical, and default) versus rectangular like a histogram.
#' (Note: as of the time shape "hist" was added, combination of jittered points is not supported by the \code{\link[ggridges]{stat_binline}} that dittoSeq relies on.)
#' @param ridgeplot.bins Integer which sets how many chunks to break the x-axis into when \code{ridgeplot.shape = "hist"}.
#' Overridden by \code{ridgeplot.binwidth} when that input is provided.
#' @param ridgeplot.binwidth Integer which sets the width of chunks to break the x-axis into when \code{ridgeplot.shape = "hist"}.
#' Takes precedence over \code{ridgeplot.bins} when provided.
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.title String or \code{NULL}, sets the title for the main legend which includes colors and data representations.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot (\code{p}) and data (\code{data}).
#' @param ... arguments passed to dittoPlot by dittoRidgePlot, dittoRidgeJitter, and dittoBoxPlot wrappers.
#' Options are all the ones above.
#'
#' @inheritParams gene
#'
#' @return a ggplot where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively when \code{data.out=TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the ggplot where additional data will be displayed when the cursor is hovered over jitter points.
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided, that data).
#' On the discrete axis, data will be grouped by the metadata given to \code{group.by} and colored by the metadata given to \code{color.by}.
#' The \code{assay} and \code{slot} inputs can be used to change what expression data is used when displaying gene expression.
#' If a set of cells to use is indicated with the \code{cells.use} input, the data is subset to include only those cells before plotting.
#'
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Inclusion of \code{"ridgeplot"} overrides \code{"boxplot"} and \code{"vlnplot"} presence and changes the plot to be horizontal.
#'
#' When \code{split.by} is provided the name of a metadata containing discrete data, separate plots will be produced representing each of the distinct groupings of the split.by data.
#'
#' \code{dittoRidgePlot}, \code{dittoRidgeJitter}, and \code{dittoBoxPlot} are included as wrappers of the basic \code{dittoPlot} function
#' that simply change the default for the \code{plots} input to be \code{"ridgeplot"}, \code{c("ridgeplot","jitter")}, or \code{c("boxplot","jitter")},
#' to make such plots even easier to produce.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' The \code{plots} argument determines the types of \strong{data representation} that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' 
#' Each plot type has specific associated options which are controlled by variables that start with their associated string.
#' For example, all jitter adjustments start with "\code{jitter.}", such as \code{jitter.size} and \code{jitter.width}.
#'
#' Inclusion of \code{"ridgeplot"} overrides \code{"boxplot"} and \code{"vlnplot"} presence and changes the plot to be horizontal.
#'
#' Additionally:
#'
#' \itemize{
#' \item \strong{Colors can be adjusted} with \code{color.panel}.
#' \item \strong{Subgroupings:} \code{color.by} can be utilized to split major \code{group.by} groupings into subgroups.
#' When this is done in y-axis plotting, dittoSeq automatically ensures the centers of all geoms will align,
#' but users will need to manually adjust \code{jitter.width} to less than 0.5/num_subgroups to avoid overlaps.
#' There are also three inputs through which one can use to control geom-center placement, but the easiest way to do all at once so is to just adjust \code{vlnplot.width}!
#' The other two: \code{boxplot.position.dodge}, and \code{jitter.position.dodge}.
#' \item \strong{Line(s) can be added} at single or multiple value(s) by providing these values to \code{add.line}.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' \item \strong{Titles and axes labels} can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The \strong{legend can be hidden} by setting \code{legend.show = FALSE}.
#' \item \strong{y-axis zoom and tick marks} can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item \strong{x-axis labels and groupings} can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned on/off with \code{x.labels.rotate = TRUE/FALSE}.
#' \item \strong{Shapes used} in conjunction with \code{shape.by} can be adjusted with \code{shape.panel}.
#' \item Single or multiple \strong{additional per-cell features can be retrieved} and stashed within the underlying data using \code{extra.vars}.
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
#' dittoPlot(object = myRNA, var = "gene1", group.by = "conditions",
#'     color.by = "timepoint")
#'
#' # Update the 'plots' input to change / reorder the data representations
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"))
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("ridgeplot", "jitter"))
#'
#' ### Provided wrappers enable certain easy adjustments of the 'plots' parameter.
#' # Quickly make a Boxplot
#' dittoBoxPlot(myRNA, "gene1", group.by = "timepoint")
#' # Quickly make a Ridgeplot, with or without jitter
#' dittoRidgePlot(myRNA, "gene1", group.by = "timepoint")
#' dittoRidgeJitter(myRNA, "gene1", group.by = "timepoint")
#'
#' ### Additional Functionality
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
#' # Multiple genes or continuous metadata can also be plotted by giving them as
#' #   a vector to 'var'. One aesthetic of the plot will then be used to display
#' #   'var'-info, and you can control which (faceting / "split", x-axis grouping
#' #   / "group", or color / "color") with 'multivar.aes':
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"))
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"),
#'     multivar.aes = "group")
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"),
#'     multivar.aes = "color")
#' 
#' # For faceting, instead of using 'split.by', the target data can alternatively
#' #   be given to 'extra.var' to have it added in the underlying dataframe, then
#' #   faceting can be added manually for extra flexibility
#' dittoPlot(myRNA, "gene1", "clustering",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     extra.var = "SNP") + facet_wrap("SNP", ncol = 1, strip.position = "left")
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
    multivar.aes = c("split", "group", "color"),
    multivar.split.dir = c("col", "row"),
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
    min = NA,
    max = NA,
    xlab = "make",
    x.labels = NULL,
    x.labels.rotate = NA,
    x.reorder = NULL,
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
    do.raster = FALSE,
    raster.dpi = 300,
    jitter.size = 1,
    jitter.width = 0.2,
    jitter.color = "black",
    jitter.shape.legend.size = NA,
    jitter.shape.legend.show = TRUE,
    jitter.position.dodge = boxplot.position.dodge,
    boxplot.width = 0.2,
    boxplot.color = "black",
    boxplot.show.outliers = NA,
    boxplot.outlier.size = 1.5,
    boxplot.fill = TRUE,
    boxplot.position.dodge = vlnplot.width,
    boxplot.lineweight = 1,
    vlnplot.lineweight = 1,
    vlnplot.width = 1,
    vlnplot.scaling = "area",
    vlnplot.quantiles = NULL,
    ridgeplot.lineweight = 1,
    ridgeplot.scale = 1.25,
    ridgeplot.ymax.expansion = NA,
    ridgeplot.shape = c("smooth", "hist"),
    ridgeplot.bins = 30,
    ridgeplot.binwidth = NULL,
    add.line = NULL,
    line.linetype = "dashed",
    line.color = "black",
    legend.show = TRUE,
    legend.title = "make",
    data.out = FALSE) {

    ridgeplot.shape <- match.arg(ridgeplot.shape)
    multivar.aes <- match.arg(multivar.aes)
    multivar.split.dir <- match.arg(multivar.split.dir)
    
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
    xlab <- .leave_default_or_null(xlab,
        default = group.by,
        null.if = multivar.aes=="group" && length(var)>1 && length(var) != length(all.cells))
    main <- .leave_default_or_null(main, var,
        null.if = !(length(var)==1 && is.character(var)))
    legend.title <- .leave_default_or_null(legend.title, var,
        null.if = is.null(shape.by))

    # Grab the data
    gather_out <- .dittoPlot_data_gather(object, var, group.by, color.by,
        c(shape.by,split.by,extra.vars), cells.use, assay, slot, adjustment,
        swap.rownames, do.hover, hover.data, x.reorder, x.labels,
        split.by, multivar.aes, multivar.split.dir)
    Target_data <- gather_out$Target_data
    split.by <- gather_out$split.by

    # Make the plot
    p <- ggplot(Target_data, aes(fill=.data$color)) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)
    if(!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, Target_data, plots, xlab, ylab, shape.by, jitter.size,
            jitter.width, jitter.color, shape.panel, jitter.shape.legend.size,
            jitter.shape.legend.show, jitter.position.dodge,
            do.raster, raster.dpi,
            boxplot.width, boxplot.color, boxplot.show.outliers,
            boxplot.outlier.size, boxplot.fill,
            boxplot.position.dodge, boxplot.lineweight,
            vlnplot.lineweight, vlnplot.width, vlnplot.scaling,
            vlnplot.quantiles,
            add.line, line.linetype, line.color,
            x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, Target_data, plots, xlab, ylab, jitter.size, jitter.color,
            jitter.shape.legend.size, jitter.shape.legend.show,
            ridgeplot.lineweight, ridgeplot.scale, ridgeplot.ymax.expansion,
            ridgeplot.shape, ridgeplot.bins, ridgeplot.binwidth, add.line,
            line.linetype, line.color, x.labels.rotate, do.hover, color.panel,
            colors, y.breaks, min, max)
    }
    # Extra tweaks
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, split.adjust)
    }
    
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    if (do.hover) {
        p <- .warn_or_jitter_plotly(p, plots)
    }
    
    # DONE. Return the plot +/- data
    if (data.out) {
        list(
            p = p,
            data = Target_data)
    } else {
        p
    }
}

#' @describeIn dittoPlot Plots continuous data for customizeable cells'/samples' groupings horizontally in a density representation
#' @export
dittoRidgePlot <- function(..., plots = c("ridgeplot")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot dittoRidgePlot, but with jitter overlaid
#' @export
dittoRidgeJitter <- function(..., plots = c("ridgeplot", "jitter")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot Plots continuous data for customizeable cells'/samples' groupings in boxplot form
#' @export
dittoBoxPlot <- function(..., plots = c("boxplot","jitter")){ dittoPlot(..., plots = plots) }

.dittoPlot_add_data_y_direction <- function(
    p, Target_data, plots, xlab, ylab, shape.by,
    jitter.size, jitter.width, jitter.color,shape.panel,
    jitter.shape.legend.size, jitter.shape.legend.show, jitter.position.dodge,
    do.raster, raster.dpi,
    boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.outlier.size,
    boxplot.fill, boxplot.position.dodge, boxplot.lineweight,
    vlnplot.lineweight, vlnplot.width, vlnplot.scaling, vlnplot.quantiles,
    add.line, line.linetype, line.color,
    x.labels.rotate, do.hover, y.breaks, min, max,
    object) {
    # This function takes in a partial dittoPlot ggplot object without any data
    # overlay, and parses adding the main data visualizations.
    # Adds plots based on what is requested in plots, ordered by their order.

    # Now that we know the plot's direction, set direction & y-axis limits
    p <- p + aes(x = .data$grouping, y = .data$var.data)
    
    if (!is.null(y.breaks)) {
        p <- p + scale_y_continuous(breaks = y.breaks)
    }
    if (!is.na(min) || !is.na(max)) {
        p <- p + coord_cartesian(ylim=c(min,max))
    }

    # Add Plots
    for (i in seq_along(plots)) {
        if (plots[i] == "vlnplot") {
            p <- p + geom_violin(
                linewidth = vlnplot.lineweight,
                width = vlnplot.width,
                scale = vlnplot.scaling,
                draw_quantiles = vlnplot.quantiles,
                na.rm = TRUE)
        }

        if (plots[i] == "boxplot") {
            boxplot.args <- list(
                width = boxplot.width,
                color = boxplot.color,
                lwd = boxplot.lineweight,
                alpha = ifelse(boxplot.fill, 1, 0),
                position = position_dodge(width = boxplot.position.dodge),
                outlier.size = boxplot.outlier.size,
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
            
            # Create geom_jitter() arguments
            jitter.args <- list(
                position = position_jitterdodge(
                      jitter.width = jitter.width,
                      jitter.height = 0,
                      dodge.width = jitter.position.dodge,
                      seed = NA
                ),
                size=jitter.size,
                color = jitter.color)
            
            geom_for_jitter <- geom_jitter
            if (do.raster) {
                .error_if_no_ggrastr()
                geom_for_jitter <- ggrastr::geom_jitter_rast
                jitter.args$raster.dpi <- raster.dpi
            }
            
            jitter.aes <- aes()
            if (do.hover) {
                jitter.aes <- modifyList(jitter.aes, aes(text = .data$hover.string))
            }
            
            #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
            if (!is.null(shape.by) && isMeta(shape.by, object)) {
                
                # Set shape in aes & set scales/theming.
                jitter.aes <- modifyList(jitter.aes, aes(shape = .data[[shape.by]]))
                
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
            
            jitter.args$mapping <- jitter.aes
            
            if (do.hover) {
                p <- p + suppressWarnings(do.call(geom_for_jitter, jitter.args))
            } else {
                p <- p + do.call(geom_for_jitter, jitter.args)
            }
        }
    }

    # Add labels and, if requested, lines
    p <- p + xlab(xlab) + ylab(ylab)
    if (is.na(x.labels.rotate) || x.labels.rotate) {
        p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1))
    }
    if (!is.null(add.line)) {
        p <- p + geom_hline(yintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

#' @importFrom ggridges geom_density_ridges2
.dittoPlot_add_data_x_direction <- function(
    p, Target_data, plots, xlab, ylab, jitter.size, jitter.color,
    jitter.shape.legend.size, jitter.shape.legend.show,
    ridgeplot.lineweight, ridgeplot.scale,
    ridgeplot.ymax.expansion, ridgeplot.shape, ridgeplot.bins,
    ridgeplot.binwidth, add.line, line.linetype, line.color,
    x.labels.rotate, do.hover, color.panel, colors, y.breaks, min, max) {
    #This function takes in a partial dittoPlot ggplot object without any data overlay, and parses adding the main data visualizations.

    # Now that we know the plot's direction, set direction & "y"-axis limits
    p <- p + aes(x = .data$var.data, y = .data$grouping)
    
    if (!is.null(y.breaks)) {
        p <- p + scale_x_continuous(breaks = y.breaks)
    }
    if (!is.na(min) || !is.na(max)) {
        p <- p + coord_cartesian(xlim=c(min,max))
    }
    
    # For stylistic issues with plotting defaults, also adjust grouping-axis limits
    if (is.na(ridgeplot.ymax.expansion)) {
        num_groups <- length(unique(Target_data$grouping))
        # From 0.6 to 0.1 between 4 to 12 groups and 0.05 by 34 groups.
        set_exp <- stats::approxfun(
            x=c(3,12,34), y=c(0.6, 0.1, 0.05), yleft = 0.6, yright = 0.05)
        ridgeplot.ymax.expansion <- set_exp(num_groups)
    }
    p <- p + scale_color_manual(values=color.panel[colors]) +
        scale_y_discrete(expand = expansion(mult=c(0, ridgeplot.ymax.expansion)))

    # Add ridgeplot and jitter data
    ridge.args <- list(linewidth = ridgeplot.lineweight, scale = ridgeplot.scale)
    if (ridgeplot.shape == "hist") {
        ridge.args$stat <- "binline"
        ridge.args$bins <- ridgeplot.bins
        ridge.args$binwidth <- ridgeplot.binwidth
    }
    if ("jitter" %in% plots) {
        ridge.args <- c(ridge.args, jittered_points = TRUE,
            point_size = jitter.size, point_color = jitter.color)
    }
    
    p <- p + do.call(ggridges::geom_density_ridges2, ridge.args)
        
    if (!is.na(jitter.shape.legend.size)) {
        p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
    }
    if (jitter.shape.legend.show==FALSE){
        p <- p + guides(shape = "none")
    }

    # Add labels and, if requested, lines
    p <- p + xlab(ylab) + ylab(xlab)
    if (!is.na(x.labels.rotate) && x.labels.rotate) {
        p <- p + theme(axis.text.y= element_text(angle=45, hjust = 1, vjust = 1))
    }
    if (!is.null(add.line)) {
        p <- p + geom_vline(xintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

.warn_or_jitter_plotly <- function(p, plots) {
    if ("ridgeplot" %in% plots) {
        warning("'do.hover = TRUE' request ignored because plotly does not support ridgeplots.")
    } else {
        .error_if_no_plotly()
        # Add hover.text to jitter, else just convert.
        if ("jitter" %in% plots) {
            p <- plotly::ggplotly(p, tooltip = "text")
        } else {
            p <- plotly::ggplotly(p)
        }
    }
    p
}

.dittoPlot_data_gather <- function(
    object, var, group.by, color.by,
    extra.vars, cells.use,
    assay, slot, adjustment,
    swap.rownames,
    do.hover, hover.data = c(var, extra.vars),
    x.reorder, x.labels,
    split.by, multivar.aes, multivar.split.dir) {

    all.cells <- .all_cells(object)
    
    # Support multiple genes/metadata
    if (length(var)>1 && length(var) != length(all.cells)) {
        Target_data <- do.call(rbind, lapply(
            var, function(this.var) {
                col <- switch(
                    multivar.aes, "split"="var", "group"="grouping", "color"="color")
                this.out <- .dittoPlot_data_gather_inner(
                    object, this.var, group.by, color.by, extra.vars, cells.use,
                    all.cells, assay, slot, adjustment, swap.rownames, do.hover,
                    hover.data, x.reorder, x.labels
                )
                this.out[[col]] <- this.var
                this.out
            }
        ))
        if (multivar.aes == "split") {
            split.by <- .multivar_adjust_split_by(
                split.by, multivar.split.dir, multivar.col.name = "var")
        }
    } else {
        # Single var
        Target_data <- .dittoPlot_data_gather_inner(
            object, var, group.by, color.by, extra.vars, cells.use,
            all.cells, assay, slot, adjustment, swap.rownames, do.hover,
            hover.data, x.reorder, x.labels
        )
    }
    list(Target_data = Target_data, split.by = split.by)
}

.dittoPlot_data_gather_inner <- function(
    object, var, group.by, color.by, extra.vars, cells.use, all.cells,
    assay, slot, adjustment, swap.rownames, do.hover, hover.data,
    x.reorder, x.labels
) {
    
    # Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    
    ### Make dataframe for storing the plotting data:
    full_data <- data.frame(
        var.data = .var_OR_get_meta_or_gene(
            var, object, assay, slot, adjustment, swap.rownames),
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

    Target_data <- full_data[all.cells %in% cells.use,]
    # Reorder / Relabel grouping data
    Target_data$grouping <-
        .rename_and_or_reorder(Target_data$grouping, x.reorder, x.labels)
    
    Target_data
}

