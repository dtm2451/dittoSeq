########## dittoBarPlot: Builds a stacked bar plot to show the composition of samples / ages / 'group.by' ##########
#' Outputs a stacked bar plot to show the percent composition of samples or other cell groupings
#' @import ggplot2
#'
#' @param var String name of a metadata that contains discrete data, or a factor or vector containing such data for all cells/samples in the target \code{object}. REQUIRED.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @param group.by String representing the name of a "metadata" to use for separating the cells/samples into discrete groups. REQUIRED.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#'
#' NOTE: When \code{cells.use} is combined with \code{scale = "percent"} left out cells are not considered in calculating percentages. Percents will always total to 1.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param scale "count" or "percent". Sets whether data should be shown as raw counts or scaled to 1 and shown as a percentage.
#' @param do.hover Logical which sets whether the ggplot output should be converted to a ggplotly object with data about individual bars displayed when you hover your cursor over them.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments. Default = \code{theme_classic()}. See \code{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options.
#' @param theme ggplot theme. Default = theme_classic()
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param ylab String, sets the continuous-axis label (=y-axis for box and violin plots, x-axis for ridgeplots).
#' @param x.labels String vector which will replaceme the x-axis grouping labels.
#' The first component of \code{x.labels} sets the name for the first x-axis grouping.
#' @param x.labels.rotate Logical which sets whether the x-axis grouping labels should be rotated.  Default = FALSE = vertical labels.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of x-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of \code{x.reorder} should be these indices, but in the order that you would like them rearranged to be.
#' @param y.breaks Numeric vector which indicates the plots major gridlines. c(break1,break2,break3,etc.)
#' Note: The maximum value of this vector will set the maximum value portrayed, even if the data extends beyond this value.
#' @param main String, sets the plot title
#' @param sub String, sets the plot subtitle
#' @param var.labels.rename String vector which renames for the identities of \code{var} groupings.
#' @param var.labels.reorder Integer vector. A sequence of numbers, from 1 to the number of distinct var labels, for rearranging the order of labels' groupings within the plot.
#'
#' Method: Make a first plot without this input.
#' Then, treating the top-most grouping as index 1, and the bottom-most as index n.
#' Values of \code{var.labels.reorder} should be these indices, but in the order that you would like them rearranged to be.
#' @param legend.show Logical which sets whether the legend should be displayed. Default = TRUE.
#' @param legend.title String which adds a title to the legend.
#' @param data.out Logical which to output a dataframe containing the underlying data instead of outputing the plot itself.
#' @return A ggplot plot where discrete data, grouped by sample, condition, cluster, etc. on the x-axis, is shown on the y-axis as either counts or percent-of-total-per-grouping in a stacked barplot.
#' Alternatively, if \code{data.out = TRUE} is added, outputs the underlying data for such a plot.
#' Alternatively, if \code{do.hover = TRUE} is added, outputs a plotly conversion of such a ggplot in which underlying data can be retrieved upon hovering the cursor over the plot.
#' @examples
#' library(Seurat)
#' pbmc <- pbmc_small
#'
#' dittoBarPlot("RNA_snn_res.0.8", object = "pbmc", group.by = "ident")
#' dittoBarPlot("RNA_snn_res.0.8", object = "pbmc", group.by = "ident",
#'     scale = "count")
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand,
#'   # the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' dittoBarPlot("RNA_snn_res.0.8", group.by = "ident")
#'
#' # Accessing underlying data:
#' # as dataframe
#' dittoBarPlot("RNA_snn_res.0.8", object = "pbmc", group.by = "ident",
#'     data.out = TRUE)
#' # through hovering the cursor over the relevant parts of the plot
#' dittoBarPlot("RNA_snn_res.0.8", object = "pbmc", group.by = "ident",
#'     do.hover = TRUE)
#' @export

dittoBarPlot <- function(
    var, object = DEFAULT, group.by = "Sample", scale = c("percent", "count"),
    cells.use = NULL, data.out = FALSE, do.hover = FALSE,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    y.breaks = NA, var.labels.rename = NULL, var.labels.reorder = NULL,
    x.labels = NULL, x.labels.rotate = TRUE, x.reorder = NULL,
    theme = theme_classic(),
    xlab = group.by, ylab = "make", main = "make", sub = NULL,
    legend.show = TRUE, legend.title = NULL) {

    #Turn the object into a "name" if a full object was given
    if (typeof(object)=="S4") { object <- deparse(substitute(object)) }

    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)
    scale = match.arg(scale)

    # Create data.frame
    y.var <- as.character(
        .var_OR_get_meta_or_gene(var, object)[all.cells %in% cells.use])
    x.var <- as.character(
        .var_OR_get_meta_or_gene(group.by, object)[all.cells %in% cells.use])

    # Decide titles
    if (!(is.null(ylab))) {
        if (ylab == "make") {
            ylab.start <- ifelse(scale=="count", "Number of ", "Percent of ")
            ylab <- paste0(
                ylab.start,
                ifelse(grepl("RNAseq",.class_of(object)), "samples", "cells"))
        }
    }
    if (!(is.null(main)) && main=="make"){
        if (length(var)==1) {
            main <- var
        } else {
            main <- NULL
        }
    }

    # Create dataframe
    data <- data.frame(
        count = as.vector(data.frame(table(y.var, x.var))))
    names(data) <- c("label", "grouping", "count")
    data$label.count.total <- rep(
        as.vector(table(x.var)),
        each = length(levels(as.factor(y.var))))
    data$percent <- data$count / data$label.count.total
    data$grouping <- .rename_and_or_reorder(data$grouping, x.reorder, x.labels)
    data$label <- .rename_and_or_reorder(
        data$label, var.labels.reorder, var.labels.rename)
    if (scale == "percent") {
        y.show <- "percent"
    } else {
        y.show <- "count"
    }

    # Add hover info
    if (do.hover) {
        hover.data <- data[,names(data) %in% c("label", "count", "percent")]
        names(hover.data)[1] <- var
        # Make hover srtings, "data.type: data" \n "data.type: data"
        data$hover.string <- .make_hover_strings_from_df(hover.data)
    }

    #Build Plot
    p <- ggplot(
        data=data,
        aes_string(x = "grouping", y= y.show, fill = "label")) +
        theme + xlab(xlab) + ylab(ylab) + ggtitle(main, subtitle = sub) +
        scale_fill_manual(name = legend.title, values = color.panel[colors]) +
        if (x.labels.rotate) {
            theme(axis.text.x= element_text(
                angle=45, hjust = 1, vjust = 1, size=12))
        }
        #Add the bars.
        if(do.hover){
          p <- p + geom_col(aes_string(text = "hover.string"))
        } else {
          p <- p + geom_col()
        }
    # Set y-axis scaling
    if (is.na(y.breaks[1]) && scale == "percent") {
        y.breaks <- c(0,0.5,1)
    }
    if (!is.na(y.breaks[1])) {
        p <- p + scale_y_continuous(
              breaks= y.breaks, limits = c(0,max(y.breaks)))
    }

    if (!legend.show) {
        p <- remove_legend(p)
    }

    #DONE. Return the plot
    if (data.out) {
        return(data)
    } else {
        if (do.hover) {
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
          return(p)
        }
    }
}
