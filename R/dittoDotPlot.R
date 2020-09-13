#' Compact plotting of per group summaries for expression of multiple features
#' 
#' @param vars String vector (example: \code{c("gene1","gene2","gene3")}) which selects which variables, typically genes, to extract from the object, summarize across groups, and add to the plot
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param summary.fxn.color,summary.fxn.size A function which sets how color or size will be used to summarize variables' data for each group.
#' Defaults:
#' Any function can be used as long as it takes in a numeric vector and returns a single numeric value.
#' The name of the function becomes the default legned title for that parameter.
#' Alternative examples: \code{\link{median}}, \code{\link{max}}, \code{function (x) sum(x!=0)/length(x)}.
#' @param adjustment When plotting gene expression (or antibody, or other forms of counts data), should that data be adjusted altogether before splitting into groups?
#' \itemize{
#' \item{"relative.to.max": DEFAULT, divided by the maximum expression value to give percent of max values between [0,1]}
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{NULL: no adjustment}
#' }
#' @param scale String which sets whether mean expression values per group should be centered and scaled. 
#' @param size Number which sets the dot size representing 100\% expression in every sample.
#' @param min.percent,max.percent Numbers between 0 and 1 which sets the minimum and maaximum percent expression to show.
#' When set to NA, the minimum/mamximum of the data are used.
#' @param ylab String which sets the y axis label.
#' Set to \code{NULL} to remove.
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param y.labels.rotate Logical which sets whether the labels should be rotated.
#' @param y.labels String vector, c("label1","label2","label3",...) which overrides the names of the samples/groups.  NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param y.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of y-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the bottom-most grouping as index 1, and the top-most as index n,
#' values of y.reorder should be these indices, but in the order that you would like them rearranged to be.
#' @param legend.color.title,legend.size.title String or \code{NULL}, sets the title displayed above legend keys.
#' Defaults to the namee of the relevant \code{summary.fxn}.
#' 
#' @inheritParams dittoPlotVarsAcrossGroups
#' @inheritParams dittoScatterPlot
#' 
#' @return a ggplot object where dots of different colors and sizes summarize continuous data for multiple features (rows) per multiple groups (columns)
#'
#' Alternatively when \code{data.out = TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the plot where additional data will be displayed when the cursor is hovered over the dots.
#' 
#' @details
#' FILLER 
#' Generally, this function will output a dittoPlot, grouped by sample, age, cluster, etc., where each data point represents the summary (typically \code{mean}), across each group, of individual variable's expression, but variables can be genes or metadata.
#'
#' The data for each element of \code{vars} is obtained.
#' When elements are genes/features, \code{assay} and \code{slot} are utilized to determine which expression data to use,
#' and \code{adjustment} determines if and how the expression data might be adjusted.
#'
#' By default, a z-score adjustment is applied to all gene/feature \code{vars}. Note that this adjustment is applied \emph{before} cells/samples subsetting.
#'
#' x-axis groupings are then determined using \code{group.by}, and data for each variable is summarized using the \code{summary.fxn}.
#'
#' Finally, data is plotted with the data representation types in \code{plots}.
#'
#' @section Plot Customization:
#' FILLER
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
#' \item x-axis labels and groupings can be changed / reordered using \code{y.labels} and \code{y.reorder}, and rotation of these labels can be turned off with \code{y.labels.rotate = FALSE}.
#' \item Line(s) can be added at single or multiple value(s) by providing these values to \code{add.line}.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' }
#'
#' @seealso
#' \code{\link{dittoPlotVarsAcrossGroups}} for a method of summarizing expression of multiple features across distinct groups that can be better (and more compact) when the identities of the individual genes are unimportant.
#' \code{\link{dittoPlot}} and \code{\link{multi_dittoPlot}} for plotting of expression and metadata vars, each as separate plots, on a per cell/sample basis.
#'
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#' 
#' # These random data aren't very exciting, but we can at least add some zeros.
#' counts(myRNA)[1:4,1:40] <- 0
#' 
#' dittoDotPlot(
#'     myRNA, c("gene1", "gene2", "gene3", "gene4"), "clustering")
#'
#' @author Daniel Bunis
#' @export
#'
dittoDotPlot <- function(
    object,
    vars,
    group.by,
    scale = TRUE,
    cells.use = NULL,
    size = 6,
    min.percent = 0.01,
    max.percent = NA,
    min.color = "grey90",
    max.color = "#C51B7D",
    min = "make",
    max = NULL,
    summary.fxn.color = function(x) {mean(x[x>0])},
    summary.fxn.size = function(x) {mean(x>0)},
    assay = .default_assay_raw(object),
    slot = .default_slot_raw(object),
    adjustment = NULL,
    do.hover = FALSE,
    main = NULL,
    sub = NULL,
    ylab = group.by,
    y.labels = NULL,
    y.labels.rotate = TRUE,
    y.reorder = NULL,
    xlab = NULL,
    theme = theme_classic(),
    legend.show = TRUE,
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.color.title = "make",
    legend.size.title = "percent\nexpression",
    data.out = FALSE) {

    cells.use <- .which_cells(cells.use, object)
    
    # Fill defaults
    legend.color.title <- .leave_default_or_null(
        legend.color.title,
        default = ifelse(scale,"relative\nexpression","average\nexpression"))
    min <- .leave_default_or_null(
        min,
        default = if (scale) {NULL} else {0})

    # Create data table summarizing vars data for each group
    data <- .data_gather_summarize_vars_by_groups(
        object, vars, group.by,
        list(summary.fxn.color, summary.fxn.size),
        c("color", "size"),
        cells.use, assay, slot, adjustment, do.hover)
    data$grouping <-
        .rename_and_or_reorder(data$grouping, y.reorder, y.labels)
    
    if (scale) {
        for (i in vars) {
            data$pre.scale <- data$color
            data$color[data$var == i] <- scale(data$color[data$var == i])
        }
    }

    # Generate Plot
    p <- .ditto_dot_plot(
        data, do.hover, main, sub, ylab, xlab, y.labels.rotate, scale,
        min.color, max.color, min, max,
        size, min.percent, max.percent, theme,
        legend.color.title, legend.color.breaks, legend.color.breaks.labels,
        legend.size.title, legend.show)
    
    # DONE. Return
    if (data.out) {
        return(list(p = p, data = data))
    } else {
        if (do.hover) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

.ditto_dot_plot <- function(
    data,
    do.hover,
    main,
    sub,
    ylab,
    xlab,
    y.labels.rotate,
    scale,
    min.color,
    max.color,
    min,
    max,
    size,
    min.percent,
    max.percent,
    theme,
    legend.color.title,
    legend.color.breaks,
    legend.color.breaks.labels,
    legend.size.title,
    legend.show) {
    
    p <- ggplot(data,
            aes_string(x = "var", y = "grouping", color = "color", size = "size")) +
        theme +
        ggtitle(main, sub) + xlab(xlab) + ylab(ylab) +
        scale_size(
            name = legend.size.title,
            limits = c(min.percent, max.percent),
            range = c(0, size)) +
        scale_color_gradient(
            name = legend.color.title,
            low= min.color, high = max.color,
            limits = c(
                ifelse(is.null(min), min(data$color), min),
                ifelse(is.null(max), max(data$color), max)),
            breaks = legend.color.breaks,
            labels = legend.color.breaks.labels)
    
    if (do.hover) {
        p <- p + suppressWarnings(
            geom_point(aes(text = "hover.string"), na.rm = TRUE))
    } else {
        p <- p + geom_point(na.rm = TRUE)
    }
    
    if (y.labels.rotate) {
        p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1))
    }
    
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    p
}

.data_gather_summarize_vars_by_groups <- function(
    object,
    vars,
    group.by,
    summary.fxns, # list of summaries to make
    names,        # vector of what to call those summaries
    cells.use,
    assay,
    slot,
    adjustment,
    do.hover) {
    
    groupings <- meta(group.by, object)[cells.use]

    ### Grab (and adjust) vars data per cell/sample
    # rows = cells/samples
    # cols = vars
    vars_data <- .multi_var_gather_raw(
        object, vars, assay, slot, adjustment, cells.use)

    ### Summarize vars data per group
    # rows = summarized vars data
    # cols = groupings
    summarize <- function(summary.fxn) {
        data.frame(
            vapply(
                unique(groupings),
                function (this.group) {
                    vapply(
                        vars,
                        function(this.var) {
                            summary.fxn(vars_data[groupings == this.group, this.var])
                        }, FUN.VALUE = numeric(1)
                    )
                }, FUN.VALUE = numeric(length(vars))),
            row.names = vars)
    }
    summary_data <- lapply(summary.fxns, summarize)

    ### Vectorize the summary data
    # rows = individual data points; each var for group1, group2, group3,...
    data <- data.frame(
        var = rep(vars, ncol(summary_data[[1]])),
        grouping = rep(unique(groupings), each = nrow(summary_data[[1]]))
    )
    
    for (i in seq_along(summary_data)) {
        data <- cbind(data, unlist(summary_data[[i]]))
    }
    names(data) <- c("var", "grouping", names)

    if (do.hover) {
        data$hover.string <- .make_hover_strings_from_df(data)
    }

    return(data)
}

.multi_var_gather_raw <- function(
    object,
    vars,
    assay,
    slot,
    adjustment,
    cells.use) {
    
    call_meta <- isMeta(vars, object, return.values = FALSE)
    meta_vars <- vars[call_meta]
    gene_vars <- isGene(vars[!call_meta], object, assay, return.values = TRUE)
    
    if (!all(vars %in% c(meta_vars, gene_vars))) {
        stop("All 'vars' must be a metadata or gene")
    }
    
    vars_data <- if (length(meta_vars)>0) {
        as.matrix(getMetas(object, names.only = FALSE)[, meta_vars])
    } else {
        NULL
    }
    
    if (length(gene_vars)>0) {
        gene_data <- t(as.matrix(.which_data(assay,slot,object)[gene_vars, ]))
        
        if (!is.null(adjustment)) {
            if (adjustment=="z-score") {
                gene_data <- as.matrix(scale(gene_data))
            }
            if (adjustment=="relative.to.max") {
                gene_data <- apply(gene_data, 2, function(x) {x/max(x)})
            }
        }
        
        vars_data <- rbind(vars_data, gene_data)
    }
    
    vars_data[cells.use, vars]
}
