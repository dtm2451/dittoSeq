#' Outputs a stacked bar plot to show the percent composition of samples, groups, clusters, or other groupings
#' @import ggplot2
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param var String name of a metadata that contains discrete data, or a factor or vector containing such data for all cells/samples in the target \code{object}.
#' @param group.by String name of a metadata to use for separating the cells/samples into discrete groups.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' 
#' Note: When \code{cells.use} is combined with \code{scale = "percent"}, left out cells are not considered in calculating percentages. Percents will always total to 1.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default.
#' @param colors Integer vector, which sets the indexes / order, of colors from color.panel to actually use.
#' (Provides an alternative to directly modifying \code{color.panel}.)
#' @param scale "count" or "percent". Sets whether data should be shown as counts versus percentage.
#' @param do.hover Logical which sets whether the ggplot output should be converted to a ggplotly object with data about individual bars displayed when you hover your cursor over them.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_classic()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param xlab String which sets the x-axis title.
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param ylab String which sets the y-axis title.
#' Default = "make" and if left as make, a title will be automatically generated.
#' @param x.labels String vector which will replace the x-axis groupings' labels.
#' Regardless of \code{x.reorder}, the first component of \code{x.labels} sets the name for the left-most x-axis grouping.
#' @param x.labels.rotate Logical which sets whether the x-axis grouping labels should be rotated.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of x-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of \code{x.reorder} should be these indices, but in the order that you would like them rearranged to be.
#' 
#' Recommendation for advanced users: If you find yourself coming back to this input too many times, an alternative solution that can be easier long-term
#' is to make the target data into a factor, and to put its levels in the desired order: \code{factor(data, levels = c("level1", "level2", ...))}.
#' \code{\link{metaLevels}} can be used to quickly get the identities that need to be part of this 'levels' input.
#' @param y.breaks Numeric vector which sets the plot's tick marks / major gridlines. c(break1,break2,break3,etc.)
#' @param min,max Scalars which control the zoom of the plot.
#' These inputs set the minimum / maximum values of the y-axis.
#' Default = set based on the limits of the data, 0 to 1 for \code{scale = "percent"}, or 0 to maximum count for 0 to 1 for \code{scale = "count"}.
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}.
#' @param split.adjust A named list which allows extra parameters to be pushed through to the faceting function call.
#' List elements should be valid inputs to the faceting functions, e.g. `list(scales = "free")`.
#' 
#' For options, when giving 1 metadata to \code{split.by}, see \code{\link[ggplot2]{facet_wrap}},
#' OR when giving 2 metadatas to \code{split.by}, see \code{\link[ggplot2]{facet_grid}}.
#' @param main String, sets the plot title
#' @param sub String, sets the plot subtitle
#' @param var.labels.rename String vector for renaming the distinct identities of \code{var} values.
#' 
#' Hint: use \code{\link{metaLevels}} or \code{unique(<var-data>)} to assess current values.
#' @param var.labels.reorder Integer vector. A sequence of numbers, from 1 to the number of distinct \code{var} value identities, for rearranging the order of labels' groupings within the plot.
#'
#' Method: Make a first plot without this input.
#' Then, treating the top-most grouping as index 1, and the bottom-most as index n.
#' Values of \code{var.labels.reorder} should be these indices, but in the order that you would like them rearranged to be.
#' @param legend.show Logical which sets whether the legend should be displayed.
#' @param legend.title String which adds a title to the legend.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot ("p") and a data.frame ("data") containing the underlying data.
#' @param retain.factor.levels Logical which controls whether factor identities of \code{var} and \code{group.by} data should be respected.
#' Set to TRUE to faithfully reflect ordering of groupings encoded in factor levels,
#' but Note that this will also force retention of groupings that could otherwise be removed via \code{cells.use}.
#'
#' @inheritParams dittoPlot
#'
#' @return A ggplot plot where discrete data, grouped by sample, condition, cluster, etc. on the x-axis, is shown on the y-axis as either counts or percent-of-total-per-grouping in a stacked barplot.
#'
#' Alternatively, if \code{data.out = TRUE}, a list containing the plot ("p") and a dataframe of the underlying data ("data").
#'
#' Alternatively, if \code{do.hover = TRUE}, a plotly conversion of the ggplot output in which underlying data can be retrieved upon hovering the cursor over the plot.
#' @details
#' The function creates a dataframe containing counts and percent makeup of \code{var} identities for each x-axis grouping (determined by the \code{group.by} input).
#' If a set of cells/samples to use is indicated with the \code{cells.use} input, only those cells/samples are used for counts and percent makeup calculations.
#' Then, a vertical bar plot is generated (\code{ggplot2::geom_col()}) showing either percent makeup if
#' \code{scale = "percent"}, which is the default, or raw counts if \code{scale = "count"}.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item Colors can be adjusted with \code{color.panel} and/or \code{colors}.
#' \item y-axis zoom and tick marks can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item Titles can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The legend can be removed by setting \code{legend.show = FALSE}.
#' \item x-axis labels and groupings can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned off with \code{x.labels.rotate = FALSE}.
#' \item y-axis \code{var}-group labels and their order can be changed / reordered using \code{var.labels} and \code{var.labels.reorder}.
#' }
#'
#' @seealso
#' \code{\link{dittoFreqPlot}} for a data representation that focuses on pre-sample frequencies of each the \code{var}-data values individually, rather than emphasizing total makeup of samples/groups.
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

dittoBarPlot <- function(
    object,
    var,
    group.by,
    scale = c("percent", "count"),
    split.by = NULL,
    cells.use = NULL,
    retain.factor.levels = FALSE,
    data.out = FALSE,
    do.hover = FALSE,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
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
    legend.show = TRUE,
    legend.title = NULL) {
    
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
    
    # Gather data
    data <- .dittoBarPlot_data_gather(
        object, var, group.by, split.by, cells.use, x.reorder, x.labels,
        var.labels.reorder, var.labels.rename, do.hover, FALSE,
        retain.factor.levels, retain.factor.levels)
    if (scale == "percent") {
        y.show <- "percent"
    } else {
        y.show <- "count"
    }

    #Build Plot
    p <- ggplot(
        data=data,
        aes(x = .data$grouping, y= .data[[y.show]], fill = .data$label)) +
        theme + xlab(xlab) + ylab(ylab) + ggtitle(main, subtitle = sub) +
        scale_fill_manual(name = legend.title, values = color.panel[colors]) +
        if (x.labels.rotate) {
            theme(axis.text.x= element_text(
                angle=45, hjust = 1, vjust = 1))
        }
        #Add the bars.
        if(do.hover){
            p <- p + suppressWarnings(geom_col(
                aes(text = .data$hover.string)))
        } else {
            p <- p + geom_col()
        }
    # Set y-axis ticks & scaling
    if (is.na(y.breaks[1]) && scale == "percent") {
        y.breaks <- c(0,0.5,1)
    }
    if (!is.na(y.breaks[1])) {
        p <- p + scale_y_continuous(breaks = y.breaks)
    }
    if (is.null(max)) {
        max <- ifelse(scale == "percent", 1, max(data$label.count.total))
    }
    p <- p + coord_cartesian(ylim=c(min,max))
    
    ### Add extra features
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, split.adjust)
    }

    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    if (do.hover) {
        .error_if_no_plotly()
        p <- plotly::ggplotly(p, tooltip = "text")
    }

    #DONE. Return the plot
    if (data.out) {
        list(
            p = p,
            data = data)
    } else {
        p
    }
}

.dittoBarPlot_data_gather <- function(
    object, var, group.by, split.by, cells.use,
    x.reorder, x.labels,
    var.labels.reorder, var.labels.rename,
    do.hover, max.normalize = FALSE,
    retain.factor.levels.var, retain.factor.levels.group,
    make.factor.var = FALSE, keep.level.order.group = FALSE
) {
    
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)
    
    # Extract x.grouping and y.labels data
    y.var <- .var_OR_get_meta_or_gene(var, object)[all.cells %in% cells.use]
    x.var <- .var_OR_get_meta_or_gene(group.by, object)[all.cells %in% cells.use]
    
    # For pre-v1.4 compatibility
    if(!retain.factor.levels.var) {
        y.var <- as.character(y.var)
    }
    if(make.factor.var) {
        y.var <- as.factor(y.var)
    }
    x.levs <- levels(as.factor(x.var))
    if(!retain.factor.levels.group) {
        x.var <- as.character(x.var)
    }
    
    # Extract or negate-away split.by data
    facet <- "filler"
    split.data <- list()
    if (!is.null(split.by)) {
        for (by in seq_along(split.by)) {
            split.data[[by]] <- meta(split.by[by], object)[all.cells %in% cells.use]
        }
        facet <- do.call(paste, split.data)
    }
    
    # Create dataframe (per split.by group)
    data <- do.call(
        rbind,
        lapply(
            unique(facet),
            function(this_facet) {
                
                # Subset data per facet
                y.var <- y.var[facet==this_facet]
                x.var <- x.var[facet==this_facet]
                
                # Create data frame
                new <- data.frame(
                    count = as.vector(data.frame(table(y.var, x.var))))
                names(new) <- c("label", "grouping", "count")
                
                new$label.count.total.per.facet <- rep(
                    as.vector(table(x.var)),
                    each = length(levels(as.factor(y.var))))
                new$percent <- new$count / new$label.count.total.per.facet
                
                # Catch 0/0
                new$percent[is.nan(new$percent)] <- 0
                
                # Add facet info
                for (by in seq_along(split.by)) {
                    new[[split.by[by]]] <- (split.data[[by]][facet==this_facet])[1]
                }
                
                new
            }
        )
    )
    
    # max.normalization per var-label
    if (max.normalize) {
        data$count.norm <- 0
        data$percent.norm <- 0
        
        for (i in unique(data$label)) {
            this_lab <- data$label == i
            data$count.norm[this_lab] <- 
                data$count[this_lab]/max(data$count[this_lab])
            data$percent.norm[this_lab] <- 
                data$percent[this_lab]/max(data$percent[this_lab])
        }
    }
    
    # Rename/reorder
    if(keep.level.order.group){
        data$grouping <- factor(data$grouping, levels = x.levs)
    }
    data$grouping <- .rename_and_or_reorder(data$grouping, x.reorder, x.labels)
    data$label <- .rename_and_or_reorder(
        data$label, var.labels.reorder, var.labels.rename)

    # Add hover info
    if (do.hover) {
        hover.data <- data[,names(data) %in% c("label", "count", "percent")]
        names(hover.data)[1] <- var
        # Make hover strings, "data.type: data" \n "data.type: data"
        data$hover.string <- .make_hover_strings_from_df(hover.data)
    }
    
    data
}
