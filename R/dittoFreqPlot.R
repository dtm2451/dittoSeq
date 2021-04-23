#' Outputs a stacked bar plot to show the percent composition of samples, groups, clusters, or other groupings
#' @import ggplot2
#'
#' @inheritParams dittoBarPlot
#' @inheritParams dittoPlot
#' @param var String name of a metadata that contains discrete data, or a factor or vector containing such data for all cells/samples in the target \code{object}.
#' @param sample.by String name of a metadata containing which samples each cell belongs to.
#' 
#' Note that when this is not provided, there will only be one data point per grouping. Warning can be expected then for all \code{plot} options except \code{"jitter"}.
#' @param vars.use A subset of the values of \code{var}-data which should be shown.
#' If left as \code{NULL}, all values will be shown.
#' Hint: use \code{\link{metaLevels}} or \code{unique(<var-data>)} to assess options.
#' @param max.normalize Logical which sets whether the data for each var-data value (each facet) should be normalized to have the same maximum value.
#' When set to \code{TRUE}, lower frequency vars will make use of just as much plot space as higher frequency vars.
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
#' \item The \strong{legend can be hidden} by setting \code{legend.show = TRUE}.
#' \item \strong{y-axis zoom and tick marks} can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item \strong{x-axis labels and groupings} can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned off with \code{x.labels.rotate = FALSE}.
#' \item \strong{Shapes used} in conjunction with \code{shape.by} can be adjusted with \code{shape.panel}.
#' }
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA1 <- myRNA
#' colnames(myRNA) <- paste0(colnames(myRNA),"_1")
#' example(importDittoBulk, echo = FALSE)
#' myRNA <- cbind(myRNA, myRNA1)
#' myRNA <- setBulk(myRNA, FALSE) 
#' myRNA$sample <- rep(1:12, each = 10)
#' myRNA$groups <- rep(c("A", "B"), each = 60)
#' myRNA$subgroups <- rep(as.character(c(1:3,1:3,1:3,1:3)), each = 10)
#' myRNA
#' 
#' # There are three necessary inputs for this function.
#' #  var = typically this will be cell types annotations or clustering
#' #  sample.by = the name of a metadata containing sample assignment of cells.
#' #  group.by = how to group group the data on the x-axis (y-axis for ridgeplots)
#' dittoFreqPlot(myRNA,
#'     var = "clustering",
#'     sample.by = "sample",
#'     group.by = "groups")
#'     
#' # color.by can also be set differently from group.by to have the effect of
#' #  adding subgroupings:
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups",
#'     sample.by = "sample",
#'     color.by = "subgroups")
#' 
#' # The var-values shown can be subset with vars.use
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups", sample.by = "sample", color.by = "subgroups",
#'     vars.use = 1:2)
#' 
#' # Lower frequency groups can be 
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups", sample.by = "sample", color.by = "subgroups",
#'     max.normalize = TRUE)
#' 
#' @author Daniel Bunis
#' @export

dittoFreqPlot <- function(
    object,
    var,
    sample.by = NULL,
    group.by,
    color.by = group.by,
    vars.use = NULL,
    scale = c("percent", "count"),
    max.normalize = FALSE,
    plots = c("boxplot","jitter"),
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
    jitter.position.dodge = boxplot.position.dodge,
    do.raster = FALSE,
    raster.dpi = 300,
    boxplot.width = 0.4,
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
    legend.title = color.by) {
    
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
    
    # Check that sample definitions are 1:1 with groupings/colorings
    samps <- meta(sample.by, object)
    .check_1value_per_group(samps, group.by, object)
    .check_1value_per_group(samps, color.by, object)
    
    # Gather data (use split.by to ensure per- color.by & sample.by calculation)
    data <- .dittoBarPlot_data_gather(
        object, var, group.by, split.by = c(sample.by, color.by),
        cells.use, x.reorder, x.labels,
        var.labels.reorder, var.labels.rename, do.hover,
        retain.factor.levels = TRUE, max.normalize)
    
    # Subset to vars.use
    if (!is.null(vars.use)) {
        data <- data[data$label %in% vars.use,]
    }
    
    # Adjust BarPlot-ready data for dittoPlot plotter expectation
    if (scale == "percent") {
        y.show <- "percent"
    } else {
        y.show <- "count"
    }
    if (max.normalize) {
        y.show <- paste0(y.show, ".norm")
        y.breaks = NULL
        ylab <- paste("Normalized", ylab)
    }
    data$var.data <- data[[y.show]]
    
    # Set y-axis ticks
    y.breaks <- if (identical(y.breaks, NA) && scale == "percent") {
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
            jitter.color, 16, NA, TRUE, jitter.position.dodge,
            do.raster, raster.dpi,
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

.check_1value_per_group <- function(groupings, check, object) {
    
    values <- meta(check, object)
    any_non_1 <- !all(vapply(
        unique(groupings),
        function (group) {
            length(unique(values[groupings == group]))==1
        }, FUN.VALUE = logical(1)
        ))
    if (any_non_1) {
        stop("Unable to interpret ",check," input. It's data does not map 1:1 per sample.")
    }
}
