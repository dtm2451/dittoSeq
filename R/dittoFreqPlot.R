#' Outputs a stacked bar plot to show the percent composition of samples, groups, clusters, or other groupings
#' @import ggplot2
#'
#' @inheritParams dittoBarPlot
#' @inheritParams dittoPlot
#' @param var String name of a metadata that contains discrete data, or a factor or vector containing such data for all cells/samples in the target \code{object}.
#' @param sample.by String name of a metadata containing which samples each cell belongs to.
#' @param nrow,ncol Integers which set the dimensions of the facet grid.
#' @return A ggplot plot where discrete data, grouped by sample, condition, cluster, etc., is shown on the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively, if \code{data.out = TRUE}, a list containing the plot ("p") and a dataframe of the underlying data ("data").
#'
#' Alternatively, if \code{do.hover = TRUE}, a plotly conversion of the ggplot output in which underlying data can be retrieved upon hovering the cursor over the plot.
#' @details
#' The function creates a dataframe containing counts and percent makeup of \code{var} identities for each x-axis grouping (determined by the \code{group.by} input).
#' If a set of cells/samples to use is indicated with the \code{cells.use} input, only those cells/samples are used for counts and percent makeup calculations.
#' Then, a vertical bar plot is generated (\code{ggplot2::geom_col()})
#'
#' Finally, data is plotted with the data representation types in \code{plots}, showing either percent of total if
#' \code{scale = "percent"}, which is the default, or raw counts if \code{scale = "count"}.
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
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' dittoBarPlot(myRNA, "clustering", group.by = "groups")
#' dittoBarPlot(myRNA, "clustering", group.by = "groups",
#'     scale = "count")
#'
#' # Reordering the x-axis groupings to have "C" (#3) come first
#' dittoBarPlot(myRNA, "clustering", group.by = "groups",
#'     x.reorder = c(3,1,2,4))
#'
#' ### Accessing underlying data:
#' # as dataframe
#' dittoBarPlot(myRNA, "clustering", group.by = "groups",
#'     data.out = TRUE)
#' # through hovering the cursor over the relevant parts of the plot
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'     dittoBarPlot(myRNA, "clustering", group.by = "groups",
#'         do.hover = TRUE)
#'     }
#'
#' ### Previous Version Compatibility
#' # Mistakenly, dittoBarPlot used to remove factor identities entirely from the
#' #  data it used. This manifests as ignorance of a user's set orderings for
#' #  their data. That is nolonger done by default, but to recreate old plots,
#' #  restoring this behavior can be achieved with 'retain.factor.levels = FALSE'
#' # Set factor level ordering for a metadata we'll give to 'group.by'
#' myRNA$groups_reverse_levels <- factor(
#'     myRNA$groups,
#'     levels = c("D", "C", "B", "A"))
#' # dittoBarPlot will now respect this level order by default. 
#' dittoBarPlot(myRNA, "clustering", group.by = "groups_reverse_levels")
#' # But that respect can be turned off...
#' dittoBarPlot(myRNA, "clustering", group.by = "groups_reverse_levels",
#'     retain.factor.levels = FALSE)
#' 
#' @author Daniel Bunis
#' @export

dittoFreqPlot <- function(
    object,
    var,
    group.by,
    sample.by,
    color.by = group.by, # <-- think on if this makes sense!
    scale = c("percent", "count"),
    plots = c("vlnplot","jitter"),
    nrow = NULL,
    ncol = NULL,
    cells.use = NULL,
    data.out = FALSE,
    do.hover = FALSE,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    y.breaks = NA,
    min = 0,
    max = NULL,
    var.labels.rename = NULL,
    var.labels.reorder = NULL,
    x.labels = NULL,
    x.labels.rotate = TRUE,
    x.reorder = NULL,
    theme = theme_classic(),
    xlab = group.by,
    ylab = "make",
    main = "make",
    sub = NULL,
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
    retain.factor.levels = TRUE) {
    
    scale = match.arg(scale)

    # Decide titles
    if (!(is.null(ylab))) {
        if (ylab == "make") {
            ylab.start <- ifelse(scale=="count", "Number of ", "Percent of ")
            ylab <- paste0(
                ylab.start,
                ifelse(
                    isBulk(object),
                    "samples", "cells"))
        }
    }
    if (!(is.null(main)) && main=="make"){
        if (length(var)==1) {
            main <- var
        } else {
            main <- NULL
        }
    }
    
    # Gather data (use split.by to ensure per- color.by & sample.by calculation)
    data <- .dittoBarPlot_data_gather(
        object, var, group.by, split.by = c(color.by, sample.by),
        cells.use, x.reorder, x.labels,
        var.labels.reorder, var.labels.rename, do.hover, retain.factor.levels)
    
    # Adjust BarPlot-ready data for dittoPlot plotter expectation
    if (scale == "percent") {
        y.show <- "percent"
    } else {
        y.show <- "count"
    }
    data$var.data <- data[[y.show]]
    
    # Set y-axis ticks & scaling
    y.breaks <- if (is.na(y.breaks[1]) && scale == "percent") {
        c(0,0.5,1)
    } else {
        NULL
    }

    #Build Plot
    p <- ggplot(
        data=data,
        aes_string(fill = color.by)) +
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
    
    # Split by 'var' to have the desired per element effect!
    p <- .add_splitting(
        p, "label", nrow, ncol, object, cells.use)
    
    ### Add extra features
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    if (do.hover) {
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
    }
    
    # DONE. Return the plot +/- data
    if (data.out) {
        return(list(p = p, data = data))
    } else {
        return(p)
    }
}
