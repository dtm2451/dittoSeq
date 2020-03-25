#### dittoPlotVarsAcrossGroups: Generates a dittoPlot where datapoints are genes/metadata summarizes per groups instead of individual values per cells/samples.
#' Generates a dittoPlot where datapoints are genes/metadata summarizes per groups instead of individual values per cells/samples.
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param vars String vector (example: \code{c("gene1","gene2","gene3")}) which selects which variables, typically genes, to extract from the object, summarize across groups, and add to the plot
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param color.by String representing the name of a metadata to use for setting fills.
#' Great for highlighting subgroups when wanted, but it defaults to \code{group.by} so this input can be skipped otherwise.
#' Affects boxplot, vlnplot, and ridgeplot fills.
#' @param summary.fxn A function which sets how variables' data will be summarized accross the groups.
#' Default is \code{\link{mean}}, which will take the average value, but any function can be used as long as it takes in a numeric vector and returns a single numeric value.
#' Alternative examles: \code{\link{median}}, \code{\link{max}, \code{function (x) sum(x!=0)}.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{colnames(object)[USE]} OR \code{object@cell.names[USE]}.
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot".
#' Order matters: c("vlnplot", "boxplot", "jitter") will put a violin plot in the back, boxplot in the middle, and then individual dots in the front.
#' See details section for more info.
#' @param assay,slot single strings or integer that set which data to use when plotting expressin data. See \code{\link{gene}} for more information about how defaults for these are filled in when not provided.
#' @param adjustment When plotting gene expression (or antibody, or other forms of counts data), should that data be used directly or should it be adjusted to be
#' \itemize{
#' \item{"z-score": DEFAULT, scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{NULL: no adjustment, the normal method for all other ditto expression plotting}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param do.hover Logical. Default = \code{FALSE}.
#' If set to \code{TRUE} the object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.
#' The hover data works best for jitter data representations, so it is recommended to have \code{"jitter"} as the last value of the \code{plots} input when running using hover.
#'
#' Note: Currently, incompatible with RidgePlots as plotly does not support the geom.
#' @param color.panel String vector which sets the colors to draw from for plot fills.
#' Default = \code{dittoColors()}.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use.
#' (Provides an alternative to directly modifying \code{color.panel}.)
#' @param main String which sets the plot title.
#' @param sub String which sets the plot subtitle.
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_classic()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param ylab String which sets the y axis label.
#' Default = a combination of then name of the summary function + \code{adjustment} + "expression".
#' Set to \code{NULL} to remove.
#' @param y.breaks Numeric vector, a set of breaks that should be used as major gridlines. c(break1,break2,break3,etc.).
#' @param min,max Scalars which control the zoom of the plot.
#' These inputs set the minimum / maximum values of the data to show.
#' Default = set based on the limits of the data in var.
#' @param xlab String which sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param x.labels String vector, c("label1","label2","label3",...) which overrides the names of the samples/groups.  NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging the order of x-axis groupings.
#'
#' Method: Make a first plot without this input.
#' Then, treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of x.reorder should be these indices, but in the order that you would like them rearranged to be.
#' @param x.labels.rotate Logical which sets whether the labels should be rotated.
#' Default: \code{TRUE} for violin and box plots, but \code{FALSE} for ridgeplots.
#' @param add.line numeric value(s) where one or multiple line should be added
#' @param line.linetype String which sets the type of line for \code{add.line}.
#' Defaults to "dashed", but any ggplot linetype will work.
#' @param line.color String that sets the color(s) of the \code{add.line} line(s)
#' @param jitter.size Scalar which sets the size of the jitter shapes.
#' @param jitter.width Scalar that sets the width/spread of the jitter in the x direction. Ignored in ridgeplots.
#' @param jitter.color String which sets the color of the jitter shapes
#' @param boxplot.width Scalar which sets the width/spread of the boxplot in the x direction
#' @param boxplot.color String which sets the color of the lines of the boxplot
#' @param boxplot.show.outliers Logical, whether outliers should by including in the boxplot.
#' Default is \code{FALSE} when there is a jitter plotted, \code{TRUE} if there is no jitter.
#' @param boxplot.fill Logical, whether the boxplot should be filled in or not.
#' Known bug: when boxplot fill is turned off, outliers do not render.
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
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot (\code{p}) and data (\code{data}).
#'
#' Note: plotly conversion is turned off in the \code{data.out = TRUE} setting, but hover.data is still calculated.
#' @return a ggplot or plotly where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively when \code{data.out=TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the plot where additional data will be displayed when the cursor is hovered over jitter points.
#' @details
#' Generally, this function will output a dittoPlot, grouped by sample, age, cluster, etc., where each data point represents the summary (typically \code{mean}), accross each group, of individual variable's expression, but variables can be genes or metadata.
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
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' ##########
#' ### Generate some random data
#' ##########
#' # Zero-inflated Expression
#' nsamples <- 60
#' exp <- rpois(1000*nsamples, 20)
#' exp[sample(c(TRUE,TRUE,FALSE),1000*nsamples, TRUE)] <- 0
#' exp <- matrix(exp, ncol=nsamples)
#' colnames(exp) <- paste0("sample", seq_len(ncol(exp)))
#' rownames(exp) <- paste0("gene", seq_len(nrow(exp)))
#' logexp <- log2(exp + 1)
#'
#' # Metadata
#' conds <- factor(rep(c("condition1", "condition2"), each=nsamples/2))
#' timept <- rep(c("d0", "d3", "d6", "d9"), each = 15)
#' genome <- rep(c(rep(TRUE,7),rep(FALSE,8)), 4)
#' grps <- sample(c("A","B","C","D"), nsamples, TRUE)
#' clusts <- as.character(1*(tsne[,1]>0&tsne[,2]>0) +
#'                        2*(tsne[,1]<0&tsne[,2]>0) +
#'                        3*(tsne[,1]>0&tsne[,2]<0) +
#'                        4*(tsne[,1]<0&tsne[,2]<0))
#'
#' # We can add these directly during import, or after.
#' myscRNA <- importDittoBulk(x = list(counts = exp, logcounts = logexp),
#'     metadata = data.frame(conditions = conds, timepoint = timept,
#'         SNP = genome, groups = grps))
#'
#' # Pick a set of genes
#' genes <- getGenes(myscRNA)[1:30]
#'
#' dittoPlotVarsAcrossGroups(
#'     myscRNA, genes, group.by = "timepoint")
#'
#' # Color can be controlled separately from grouping with 'color.by'
#' #   Just note: all groupings must map to a single color.
#' dittoPlotVarsAcrossGroups(myscRNA, genes, "timepoint",
#'     color.by = "conditions")
#'
#' # To change it to have the violin plot in the back, a jitter on
#' #  top of that, and a white boxplot with no fill in front:
#' dittoPlotVarsAcrossGroups(myscRNA, genes, "timepoint", "conditions",
#'     plots = c("vlnplot","jitter","boxplot"),
#'     boxplot.color = "white", boxplot.fill = FALSE)
#'
#' ## Data can be summaryized in other ways by changing the summary.fxn input.
#' #  Often, it makes sense to turn off the z-score adjustment in such cases.
#' #  median
#' dittoPlotVarsAcrossGroups(myscRNA, genes, "timepoint", "conditions",
#'     summary.fxn = median,
#'     adjustment = NULL)
#' #  Percent non-zero expression
#' percent <- function(x) {sum(x!=0)/length(x)}
#' dittoPlotVarsAcrossGroups(myscRNA, genes, "timepoint", "conditions",
#'     summary.fxn = percent,
#'     adjustment = NULL)
#'
#' # To investigate the identities of outlier genes, we can turn on hovering
#' # (if the plotly package is available)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'     dittoPlotVarsAcrossGroups(
#'         myscRNA, genes, "timepoint", "conditions",
#'         do.hover = TRUE)
#' }
#'
#' @author Daniel Bunis
#' @export

dittoPlotVarsAcrossGroups <- function(
    object, vars, group.by, color.by=group.by, summary.fxn = mean,
    cells.use = NULL, plots = c("vlnplot","jitter"),
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = "z-score",
    do.hover = FALSE, main = NULL, sub = NULL,
    ylab = "make", y.breaks = NULL, min = NULL, max = NULL, xlab = group.by,
    x.labels = NULL, x.labels.rotate = NA, x.reorder = NULL,
    color.panel = dittoColors(), colors = c(seq_along(color.panel)),
    theme = theme_classic(),
    jitter.size=1, jitter.width=0.2, jitter.color = "black",
    boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = NA,
    boxplot.fill =TRUE,
    vlnplot.lineweight = 1, vlnplot.width = 1, vlnplot.scaling = "area",
    ridgeplot.lineweight = 1, ridgeplot.scale = 1.25,
    add.line=NULL, line.linetype = "dashed", line.color = "black",
    legend.show = TRUE, legend.title = NULL, data.out = FALSE){

    #Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    #Establish the full list of cell/sample names
    all.cells <- .all_cells(object)

    #### Create data table
    # Summarizes data and creates vars x groupings table
    data <- .dittoPlotVarsAcrossGroups_data_gather(
        object, vars, group.by, color.by, summary.fxn, cells.use, assay, slot,
        adjustment, do.hover)
    data$grouping <-
        .rename_and_or_reorder(as.character(data$grouping),x.reorder,x.labels)
    all.genes <- ifelse(sum(!isGene(vars, object, assay))==0, TRUE, FALSE)
    ylab <- .leave_default_or_null(ylab,
        default = paste(
            deparse(substitute(summary.fxn)),
            adjustment,
            if (all.genes) {
                "expression"
            }, sep = " "))

    #####Start making the plot
    p <- ggplot(
        data,
        aes_string(x = "grouping", y = "var.data", fill = "color")) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)

    #Add data and x/y adjustments & labels
    if(!("ridgeplot" %in% plots)) {
        p <- .dittoPlot_add_data_y_direction(
            p, data, plots, xlab, ylab, NULL, jitter.size, jitter.width,
            jitter.color, 16, NA, TRUE, boxplot.width, boxplot.color,
            boxplot.show.outliers, boxplot.fill, vlnplot.lineweight,
            vlnplot.width, vlnplot.scaling, add.line, line.linetype,
            line.color, x.labels.rotate, do.hover, y.breaks, min, max, object)
    } else {
        p <- .dittoPlot_add_data_x_direction(
            p, data, plots, xlab, ylab, jitter.size, jitter.color,
            NA, TRUE, ridgeplot.lineweight, ridgeplot.scale, add.line,
            line.linetype, line.color, x.labels.rotate, do.hover, color.panel,
            colors, y.breaks, min, max)
    }
    #Remove legend, if warrented
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    #DONE. Return the plot or data
    if (data.out) {
        return(list(p = p, data = data))
    } else {
        if (do.hover & ("jitter" %in% plots)) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

.dittoPlotVarsAcrossGroups_data_gather <- function(
    object, vars, group.by = "Sample", color.by = group.by,
    summary.fxn = mean, cells.use = NULL, assay, slot, adjustment,
    do.hover = FALSE) {

    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    # Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    # Establish the full list of cell/sample names
    all.cells <- .all_cells(object)

    ### Grab vars and grouping/color data as a dataframes
    vars_data <- data.frame(
        vapply(
            vars,
            function(this)
                .var_OR_get_meta_or_gene(this,object,assay,slot,adjustment),
            FUN.VALUE = numeric(length(all.cells))),
        row.names = all.cells)[cells.use,]
    names(vars_data) <- vars

    groupings <- as.character(meta(group.by, object)[all.cells %in% cells.use])
    colors.data <- as.character(meta(color.by, object)[all.cells %in% cells.use])
    #### Ensure that there are no ambiguities between group.by and color.by
    if (color.by != group.by) {
        error <- sum(!vapply(
            metaLevels(group.by, object, cells.use),
            function (group) {
                length(levels(as.factor(colors.data[groupings == group])))==1
            }, FUN.VALUE = logical(1))
            )!=0
        if (error) {
            stop("Unable to interpret 'color.by' input. 'group.by' sets must map within the same 'color.by' sets.")
        }
    }

    ### Make summary data per var by cell/sample groupings
    summary_data <- data.frame(
        vapply(
            levels(as.factor(groupings)),
            function (this.group) {
                vapply(
                    vars,
                    function(this.var) {
                        summary.fxn(vars_data[groupings == this.group, this.var])
                    }, FUN.VALUE = numeric(1)
                )
            }, FUN.VALUE = numeric(length(vars))),
        row.names = vars
    )

    ### Vectorize the summary_data and add var.names and grouping/coloring
    data <- data.frame(
        var = rep(vars, ncol(summary_data)),
        var.data = unlist(summary_data),
        grouping = rep(
            levels(as.factor(groupings)), each = nrow(summary_data))
    )
    data$color <- colors.data[match(data$grouping, groupings)]

    if (do.hover) {
        hover.data <- data
        names(hover.data)[2] <- "value"
        data$hover.string <- .make_hover_strings_from_df(data)
    }

    return(data)
}
