#' Plot cell type/cluster/identity frequencies per sample and per grouping
#'
#' @param var String name of a metadata that contains discrete data, or a factor or vector containing such data for all cells/samples in the target \code{object}.
#' @param sample.by String name of a metadata containing which samples each cell belongs to.
#' 
#' Note that when this is not provided, there will only be one data point per grouping.
#' A warning can be expected then for all \code{plots} options except \code{"jitter"},
#' but this can be a useful exercise when simply trying to quantify cell type frequency fluctuations among 2 particular metadata (given to group.by & color.by) combinations.
#' @param vars.use String or string vector naming a subset of the values of \code{var}-data which should be shown.
#' If left as \code{NULL}, all values are shown.
#' 
#' Hint: use \code{\link{metaLevels}} or \code{unique(<var-data>)} to assess options.
#' 
#' Note: When \code{var.labels.rename} is jointly utilized to update how the \code{var}-values are shown, the \strong{updated} values must be used.
#' @param var.labels.rename String vector for renaming the distinct identities of \code{var}-values.
#' 
#' Hint: use \code{\link{metaLevels}} or \code{unique(<var-data>)} to assess current values.
#' @param var.labels.reorder Integer vector. A sequence of numbers, from 1 to the number of distinct \code{var}-value identities, for rearranging the order of facets within the plot space.
#'
#' Method: Make a first plot without this input.
#' Then, treating the top-left-most grouping as index 1, and the bottom-right-most as index n.
#' Values of \code{var.labels.reorder} should be these indices, but in the order that you would like them rearranged to be.
#' @param max.normalize Logical which sets whether the data for each \code{var}-data value (each facet) should be normalized to have the same maximum value.
#' 
#' When set to \code{TRUE}, lower frequency \code{var}-values will make use of just as much plot space as higher frequency vars.
#' @param split.nrow,split.ncol Integers which set the dimensions of the facet grid. (the \code{split.nrow} and \code{split.ncol} equivalent of other functions)
#' @param split.adjust A named list which allows extra parameters to be pushed through to the faceting function call.
#' List elements should be valid inputs to the faceting functions, e.g. `list(scales = "free")`.
#' 
#' Faceting for this dittoFreqPlot is always by the \code{var}-data, so see \code{\link[ggplot2]{facet_wrap}} for options.
#' @param ylab String, sets the continuous-axis label (=y-axis for box and violin plots, x-axis for ridgeplots).
#' Default = "make" and if left as make, a title will be automatically generated.
#'
#' @inheritParams dittoPlot
#' @inheritParams dittoBarPlot
#'
#' @return A ggplot plot where frequencies of discrete data, grouped by sample, condition, etc., is shown on the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively, if \code{data.out = TRUE}, a list containing the plot ("p") and a dataframe of the underlying data ("data").
#'
#' Alternatively, if \code{do.hover = TRUE}, a plotly conversion of the ggplot output in which underlying data can be retrieved upon hovering the cursor over the plot.
#' @details
#' The function creates a dataframe containing counts and percent makeup of \code{var} identities per sample if \code{sample.by} is given, or per group if only \code{group.by} is given.
#' \code{color.by} can optionally be used to add subgroupings to calculations and ultimate plots, or to convey super-groups of \code{group.by} groupings. 
#' 
#' Typically, \code{var} will be pointed to clustering or cell type annotations, but in truth it can be given any discrete data.
#' 
#' If a set of cells to use is indicated with the \code{cells.use} input, only those cells/samples are used for counts and percent makeup calculations.
#' 
#' If a set of \code{var}-values to show is indicated with the \code{vars.use} input, the data.frame is trimmed at the end to include only corresponding rows.
#' 
#' If \code{max.normalized} is set to \code{TRUE}, counts and percent data are transformed to a 0-1 scale, which makes better use of white space for lower frequency \code{var}-values.
#'
#' Either percent of total (\code{scale = "percent"}), which is the default, or counts (if \code{scale = "count"})
#' data is then (gg)plotted with the data representation types in \code{plots} by utilizing the same machinery as \code{\link{dittoPlot}}.
#' Faceting by \code{var}-data values is utilized to achieve per \code{var}-value (e.g. cluster or cell type) granularity.
#' 
#' See below for additional customization options!
#' 
#' @section Calculation Details:
#' The function is restricted in that each samples' cells, indicated by the unique values of \code{sample.by}-data, must exist within single \code{group.by} and \code{color.by} groupings.
#' Thus, in order to ensure all valid \code{var}-data composition data points are generated, prior to calculations... \itemize{
#' \item \code{var}-data are ensured to be a factor, which ensures a calculation will be run for every \code{var}-value (a.k.a. cell type or cluster)
#' \item \code{group.by}-data and \code{color-by}-data are treated as non-factor data, which ensures that calculations are run only for the groupings that each sample is associated with.
#' }
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
#' \item The \strong{legend can be hidden} by setting \code{legend.show = FALSE}.
#' \item \strong{y-axis zoom and tick marks} can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item \strong{x-axis labels and groupings} can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned on/off with \code{x.labels.rotate = TRUE/FALSE}.
#' }
#' 
#' @seealso
#' \code{\link{dittoBarPlot}} for a data representation that emphasizes total makeup of samples/groups rather than focusing on the \code{var}-data values individually.
#'
#' @examples
#' # Establish some workable example data
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
#' # There are three main inputs for this function, in addition to 'object'.
#' #  var = typically this will be cell types annotations or clustering
#' #  sample.by = the name of a metadata containing sample assignment of cells.
#' #  group.by = how to group the data on the x-axis (y-axis for ridgeplots)
#' dittoFreqPlot(myRNA,
#'     var = "clustering",
#'     sample.by = "sample",
#'     group.by = "groups")
#'     
#' # 'color.by' can also be set differently from 'group.by' to have the effect
#' #  of highlighting supersets or subgroupings:
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups",
#'     sample.by = "sample",
#'     color.by = "subgroups")
#'
#' # The var-values shown can be subset with 'vars.use'
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups", sample.by = "sample", color.by = "subgroups",
#'     vars.use = 1:2)
#' 
#' # Lower frequency groups can be expanded to use the entire y-axis by:
#' #  turning on 'max.normalize'-ation: 
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups", sample.by = "sample", color.by = "subgroups",
#'     max.normalize = TRUE)
#' #  or by setting y-scale limits to be set by the contents of facets:
#' dittoFreqPlot(myRNA, "clustering",
#'     group.by = "groups", sample.by = "sample", color.by = "subgroups",
#'     split.adjust = list(scales = "free_y")) 
#' 
#' # Data representations can also be selected and reordered with the 'plots'
#' #  input, and further adjusted with inputs applying to each representation.
#' dittoFreqPlot(myRNA,
#'     var = "clustering", sample.by = "sample", group.by = "groups",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     vlnplot.lineweight = 0.2,
#'     boxplot.fill = FALSE,
#'     boxplot.lineweight = 0.2)
#' 
#' # Finally, 'sample.by' is not technically required. When not given, a
#' #  single-datapoint of overall composition stats will be shown for each
#' #  grouping.
#' #  Just note, all data representation other than "jitter" will complain
#' #  due to there only being the one datapoint per group.
#' dittoFreqPlot(myRNA,
#'     var = "clustering", group.by = "groups", color.by = "subgroups",
#'     plots = "jitter") 
#' 
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
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
    cells.use = NULL,
    data.out = FALSE,
    do.hover = FALSE,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    y.breaks = NULL,
    min = 0,
    max = NA,
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
    legend.title = color.by) {
    
    scale = match.arg(scale)
    ridgeplot.shape <- match.arg(ridgeplot.shape)

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
    if (!is.null(sample.by)) {
        samps <- meta(sample.by, object)
        .check_1value_per_group(samps, group.by, "group.by", object)
        .check_1value_per_group(samps, color.by, "color.by", object)
    }
    
    # Gather data (use split.by to ensure per- color.by & sample.by calculation)
    data <- .dittoBarPlot_data_gather(
        object, var, group.by, split.by = c(sample.by, color.by),
        cells.use, x.reorder, x.labels,
        var.labels.reorder, var.labels.rename, do.hover, max.normalize,
        TRUE, FALSE, TRUE, TRUE)
    
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

    #Build Plot
    p <- ggplot(
        data=data,
        aes(fill = .data[[color.by]])) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)

    # Add data to plot
    if (!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, data, plots, xlab, ylab, NULL, jitter.size, jitter.width,
            jitter.color, 16, NA, TRUE, jitter.position.dodge,
            do.raster, raster.dpi,
            boxplot.width, boxplot.color, boxplot.show.outliers,
            boxplot.outlier.size, boxplot.fill,
            boxplot.position.dodge, boxplot.lineweight, vlnplot.lineweight,
            vlnplot.width, vlnplot.scaling, vlnplot.quantiles,
            add.line, line.linetype, line.color,
            x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, data, plots, xlab, ylab, jitter.size, jitter.color, NA, TRUE,
            ridgeplot.lineweight, ridgeplot.scale, ridgeplot.ymax.expansion,
            ridgeplot.shape, ridgeplot.bins, ridgeplot.binwidth,
            add.line, line.linetype, line.color,
            x.labels.rotate, do.hover, color.panel,
            colors, y.breaks, min, max)
    }
    
    # Split by 'var' to have the desired per element effect!
    p <- .add_splitting(
        p, "label", split.nrow, split.ncol, split.adjust)
    
    ### Add extra features
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    if (do.hover) {
        p <- .warn_or_jitter_plotly(p, plots)
    }
    
    # DONE. Return the plot +/- data
    if (data.out) {
        return(list(p = p, data = data))
    } else {
        return(p)
    }
}

.check_1value_per_group <- function(groupings, check, input.name, object) {
    
    values <- meta(check, object)
    any_non_1 <- !all(vapply(
        unique(groupings),
        function (group) {
            length(unique(values[groupings == group]))==1
        }, FUN.VALUE = logical(1)
        ))
    if (any_non_1) {
        stop("Unable to interpret '", input.name,"' with 'samples.by'. '",
             check, "' data does not map 1:1 per sample.")
    }
}
