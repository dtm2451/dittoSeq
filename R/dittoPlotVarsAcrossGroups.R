#' Generates a dittoPlot where data points are genes/metadata summaries, per groups, instead of individual values per cells/samples.
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param vars String vector (example: \code{c("gene1","gene2","gene3")}) which selects which variables, typically genes, to extract from the object, summarize across groups, and add to the plot
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param color.by String representing the name of a metadata to use for setting fills.
#' Great for highlighting subgroups when wanted, but it defaults to \code{group.by} so this input can be skipped otherwise.
#' Affects boxplot, vlnplot, and ridgeplot fills.
#' @param summary.fxn A function which sets how variables' data will be summarized across the groups.
#' Default is \code{\link{mean}}, which will take the average value, but any function can be used as long as it takes in a numeric vector and returns a single numeric value.
#' Alternative examples: \code{\link{median}}, \code{\link{max}}, or \code{function(x) mean(x!=0)}.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot".
#' Order matters: c("vlnplot", "boxplot", "jitter") will put a violin plot in the back, boxplot in the middle, and then individual dots in the front.
#' See details section for more info.
#' @param assay,slot single strings or integer that set which data to use when plotting expressin data. See \code{\link{gene}} for more information about how defaults for these are filled in when not provided.
#' @param adjustment When plotting gene expression (or antibody, or other forms of counts data), should that data be used directly or should it be adjusted to be
#' \itemize{
#' \item{"z-score": DEFAULT, centered and scaled to produce a relative-to-mean z-score representation}
#' \item{NULL: no adjustment, the normal method for all other ditto expression plotting}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param do.hover Logical. Default = \code{FALSE}.
#' If set to \code{TRUE} (and if there is a "jitter" in \code{plots}): the object will be converted to a plotly object in which underlying data about individual points will be displayed when you hover your cursor over them.
#' 
#' Note: Currently, plotly is incompatible with ridge plots as plotly does not support the geom_density_ridges2 geom.
#' @param color.panel String vector which sets the colors to draw from for plot fills.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use.
#' (Provides an alternative to directly modifying \code{color.panel}.)
#' @param main String which sets the plot title.
#' @param sub String which sets the plot subtitle.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_classic()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param ylab String which sets the y axis label.
#' Default = a combination of the name of the summary function + \code{adjustment} + "expression".
#' Set to \code{NULL} to remove.
#' @param y.breaks Numeric vector, a set of breaks that should be used as major grid lines. c(break1,break2,break3,etc.).
#' @param min,max Scalars which control the zoom of the plot.
#' These inputs set the minimum / maximum values of the data to show.
#' 
#' @inheritParams dittoPlot
#' 
#' @return a ggplot object
#'
#' Alternatively when \code{data.out = TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the plot where additional data will be displayed when the cursor is hovered over jitter points.
#' @details
#' Generally, this function will output a dittoPlot where each data point represents a gene (or metadata) rather than a cell/sample.
#' Values are the summary (\code{mean} by default) of the values for each gene or metadata requested with \code{vars}, within each group set by \code{group.by}.
#'
#' To start with, the data for each element of \code{vars} is obtained.
#' When elements are genes/features, \code{assay} and \code{slot} are utilized to determine which expression data to use,
#' and \code{adjustment} determines if and how the expression data might be adjusted.
#' By default, a z-score adjustment is applied to all gene/feature \code{vars}.
#' Note that this adjustment is applied \emph{before} cells/samples subsetting.
#'
#' x-axis groupings are then determined using \code{group.by}, and data for each variable is summarized using the \code{summary.fxn}.
#'
#' Finally, data is plotted with the data representation types in \code{plots}.
#'
#' @section Plot Customization:
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Each plot type has specific associated options which are controlled by variables that start with their associated string, ex: \code{jitter.size}.
#'
#' Inclusion of \code{"ridgeplot"} overrides boxplot and violin plot and changes the plot to be horizontal.
#'
#' \itemize{
#' \item Colors can be adjusted with \code{color.panel}.
#' \item Shapes used in conjunction with \code{shape.by} can be adjusted with \code{shape.panel}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The legend can be hidden by setting \code{legend.show = TRUE}.
#' \item y-axis zoom and tick marks can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item x-axis labels and groupings can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned off with \code{x.labels.rotate = FALSE}.
#' \item Line(s) can be added at single or multiple value(s) by providing these values to \code{add.line}.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' }
#'
#' @seealso
#' \code{\link{dittoPlot}} and \code{\link{multi_dittoPlot}} for plotting of single or mutliple expression and metadata vars, each as separate plots, on a per cell/sample basis.
#' 
#' \code{\link{dittoDotPlot}} for an alternative representation of per-group summaries of multiple vars where all vars are displayed separately, but still in a single plot.
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#'
#' # Pick a set of genes
#' genes <- getGenes(myRNA)[1:30]
#'
#' dittoPlotVarsAcrossGroups(
#'     myRNA, genes, group.by = "timepoint")
#'
#' # Color can be controlled separately from grouping with 'color.by'
#' #   Just note: all groupings must map to a single color.
#' dittoPlotVarsAcrossGroups(myRNA, genes, "timepoint",
#'     color.by = "conditions")
#'
#' # To change it to have the violin plot in the back, a jitter on
#' #  top of that, and a white boxplot with no fill in front:
#' dittoPlotVarsAcrossGroups(myRNA, genes, "timepoint",
#'     plots = c("vlnplot","jitter","boxplot"),
#'     boxplot.color = "white",
#'     boxplot.fill = FALSE)
#'
#' ## Data can be summarized in other ways by changing the summary.fxn input.
#' #  median
#' dittoPlotVarsAcrossGroups(myRNA, genes, "timepoint",
#'     summary.fxn = median,
#'     adjustment = NULL)
#' #  Percent non-zero expression ( = boring for this fake data)
#' percent <- function(x) {sum(x!=0)/length(x)}
#' dittoPlotVarsAcrossGroups(myRNA, genes, "timepoint",
#'     summary.fxn = percent,
#'     adjustment = NULL)
#'
#' # To investigate the identities of outlier genes, we can turn on hovering
#' # (if the plotly package is available)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'     dittoPlotVarsAcrossGroups(
#'         myRNA, genes, "timepoint",
#'         do.hover = TRUE)
#' }
#'
#' @author Daniel Bunis
#' @export

dittoPlotVarsAcrossGroups <- function(
    object,
    vars,
    group.by,
    color.by = group.by,
    summary.fxn = mean,
    cells.use = NULL,
    plots = c("vlnplot","jitter"),
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = "z-score",
    swap.rownames = NULL,
    do.hover = FALSE,
    main = NULL,
    sub = NULL,
    ylab = "make",
    y.breaks = NULL,
    min = NULL,
    max = NULL,
    xlab = group.by,
    x.labels = NULL,
    x.labels.rotate = NA,
    x.reorder = NULL,
    color.panel = dittoColors(),
    colors = c(seq_along(color.panel)),
    theme = theme_classic(),
    jitter.size = 1,
    jitter.width = 0.2,
    jitter.color = "black",
    do.raster = FALSE,
    raster.dpi = 300,
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
    legend.title = NULL,
    data.out = FALSE) {

    cells.use <- .which_cells(cells.use, object)
    
    ylab <- .leave_default_or_null(ylab,
        default = paste(
            deparse(substitute(summary.fxn)),
            adjustment,
            if (all(isGene(vars, object, assay))) {
                "expression"
            }, sep = " "))

    # Create data table summarizing vars data for each group
    data <- .data_gather_summarize_vars_by_groups(
        object, vars, group.by, list(summary.fxn), "var.data", cells.use,
        assay, slot, adjustment, swap.rownames, do.hover)
    
    data$grouping <-
        .rename_and_or_reorder(data$grouping, x.reorder, x.labels)
    
    data$color <- if (color.by != group.by) {
        colors.data <- meta(color.by, object)[cells.use]
        groupings <- meta(group.by, object)[cells.use]
        .check_1color_per_group(groupings, colors.data)
        
        colors.data[match(data$grouping, groupings)]
    } else {
        data$grouping
    }

    # Start making the plot
    p <- ggplot(data,
            aes_string(x = "grouping", y = "value", fill = "color")) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)

    # Add data to plot
    if (!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, data, plots, xlab, ylab, NULL, jitter.size, jitter.width,
            jitter.color, 16, NA, TRUE, do.raster, raster.dpi,
            boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
            boxplot.position.dodge, vlnplot.lineweight,
            vlnplot.width, vlnplot.scaling, add.line, line.linetype,
            line.color, x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, data, plots, xlab, ylab, jitter.size, jitter.color,
            NA, TRUE, ridgeplot.lineweight, ridgeplot.scale,
            ridgeplot.ymax.expansion, add.line, line.linetype, line.color,
            x.labels.rotate, do.hover, color.panel,
            colors, y.breaks, min, max)
    }
    
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    # DONE. Return
    if (data.out) {
        return(list(p = p, data = data))
    } else {
        if (do.hover && ("jitter" %in% plots)) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

.check_1color_per_group <- function(groupings, colors.data) {
    any_non_1 <- !all(vapply(
        unique(groupings),
        function (group) {
            length(unique(colors.data[groupings == group]))==1
        }, FUN.VALUE = logical(1)
        ))
    if (any_non_1) {
        stop("Unable to interpret 'color.by' input. Each 'group.by' set must map to a single 'color.by' set.")
    }
}
