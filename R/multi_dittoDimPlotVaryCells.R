##################### multi_dittoDimPlotVaryCells #######################
#' Generates multiple dittoDimPlots, each showing different cells, arranged into a grid.
#'
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample.
#' Alternatively, can be a vector of same length as there are cells/samples in the \code{object}.
#' Discrete or continuous data both work. REQUIRED.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, or the name of the object in "quotes". REQUIRED, unless \code{DEFAULT <- "object"} has been run.
#' @param vary.cells.meta String name of a metadata that should be used for selecting which cells to show in each "varycells" plot. REQUIRED.
#' @param vary.cells.levels The values/groupings of the \code{vary.cells.meta} metadata that should get a plot.
#' Defaults to all levels of the metadata.
#' @param show.allcells.plot Logical which sets whether an additional plot showing all of the cells should be added.
#' @param show.legend.plots Logical which sets whether or not legends should be plotted in varycells plot. Default = FALSE.
#' @param show.legend.allcells.plot Logical which sets whether or a legend should be plotted in the allcells plot. Default = FALSE.
#' @param show.legend.single Logical which sets whether to add a single legend as an additional plot. Default = TRUE.
#' @param show.titles Logical which sets whether titles should be added to the individual varycells plots
#' @param ncol,nrow Integers which set dimensions of the plot grid.
#' @param allcells.main String which adjusts the title of the allcells plot. Default = "All Cells".  Set to \code{NULL} or \code{""} to remove.
#' @param color.panel a list of colors to be used for when plotting a discrete var.
#' @param colors indexes / order of colors from color.panel to use. USAGE= changing the order of how colors are linked to specific groups
#' @param min,max Numeric values which set the values associated with the minimum and maximum colors in the legend when data represented by \code{var} is numeric.
#' @param data.type For when plotting expression data, sets the data-type slot that will be obtained. See \code{\link{gene}} for options and details. DEFAULT = "normalized".
#' @param ... additional parameters passed todittoDimPlot.
#' A few suggestions: \code{reduction.use} for setting which dimensionality reduction space to use.
#' \code{xlab} and \code{ylab} can be set to \code{NULL} to remove the axes labels and provide extra room for the data.
#' \code{size} can be used to adjust the size of the dots.
#' @param OUT.List Logical which controls whether the list of plots should be returned as a list instead of as a single grid arrangement of the plots.
#' @return multiple dittoDimPlot \code{\link[ggplot2]{ggplot}}s either arranged in a grid OR as a list
#' @details This function generates separate dittoDimPlots that show the same target data, but for distinct cells.
#' Which cells fall into which plot is controlled with the \code{vary.cells.meta} parameter.
#' When the quoted name of a metadata containing discrete groupings is given to \code{vary.cells.meta},
#' the function makes separate plots containing all cells/samples of each grouping.
#'
#' If plots for only certain groupings of cells are wanted, names of the wanted groupings can be supplied to the \code{vary.cells.levels} input.
#'
#' The function then appends a plot containing all groupings, titled as "All Cells" (unless otherwise changed with the \code{allcells.main} parameter),
#' as well as a single legend.  Either of these can be turned off with the \code{show.allcells.plot} and \code{show.legend.single} parameters.
#'
#' Plots are either  output in a grid (default) with \code{ncol} columns and \code{nrow} rows,
#' or alternatively as a simple list of ggplots if \code{OUT.List} is set to \code{TRUE}.
#' In the list, the varycells plots will be named by the \code{cells.meta.level} that they contain,
#' the allcells plot will be named "allcells" and the single legend will be named "legend".
#'
#' Either continuous or discrete \code{var} data can be displayed.
#' \itemize{
#' \item For continuous data, the range of potential values is calculated at the start, and set, so that colors represent the same values accross all plots.
#' \item For discrete data, colors used in each plot are adjusted so that colors represent the same groupings accross all plots.
#' }
#'
#' @seealso
#' \code{\link{dittoDimPlot}} for the base DimPlot plotting function
#'
#' \code{\link{multi_dittoDimPlot}} for plotting distinct \code{var}s accross plots instead of disctinct cells
#'
#' @examples
#' pbmc <- Seurat::pbmc_small
#'
#' multi_dittoDimPlotVaryCells("CD14", object = "pbmc", vary.cells.meta = "ident")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#'
#' # This function can be used to quickly scan for differences in expression
#' #   within or accross clusters/cell types by providing a gene to 'var'
#' multi_dittoDimPlotVaryCells("CD14", vary.cells.meta = "ident")
#'
#' # This function is also great for generating separate plots of each individual
#' #   element of a tsne/PCplot/umap. This can be useful to check for dispersion
#' #   of groups that might otherwise be hidden behind other cells/samples.
#' #   To do so, set 'var' and 'vary.cells.meta' the same.
#' multi_dittoDimPlotVaryCells("ident", vary.cells.meta = "ident")
#'
#' # The function can also be used to quickly visualize how separate clustering
#' #   resolutions match up to each other, or perhaps how certain conditions of
#' #   cells disperse accross clusters.
#' multi_dittoDimPlotVaryCells("RNA_snn_res.0.8", vary.cells.meta = "ident")
#'
#' # For an alternative method of viewing, and easily quantifying, how ceratin
#' #   conditions of cells disperse accross clusters, see '?dittoBarPlot'
#'
#' @author Daniel Bunis
#' @export

multi_dittoDimPlotVaryCells <- function(
    var, object = DEFAULT, vary.cells.meta,
    vary.cells.levels = meta.levels(vary.cells.meta, object),
    data.type = "normalized", min = NULL, max = NULL,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    show.titles=TRUE, show.allcells.plot = TRUE, allcells.main = "All Cells",
    show.legend.single = TRUE, show.legend.plots = FALSE,
    show.legend.allcells.plot = FALSE, nrow = NULL, ncol = NULL,
    OUT.List = FALSE, ...)
{
    if (is.character(object)) {
        object <- eval(expr = parse(text = object))
    }
    color.panel <- color.panel[colors]
    cells.meta <- meta(vary.cells.meta,object)
    #Determine if var is continuous vs discrete
    numeric <- ifelse(
        is.numeric(.var_OR_get_meta_or_gene(var[1], object, data.type)),
        TRUE,
        FALSE
    )
    if (numeric) {
        the.range <- range(.var_OR_get_meta_or_gene(var, object, data.type))
        min <- ifelse(is.null(min), the.range[1], min)
        max <- ifelse(is.null(max), the.range[2], max)
    }

    ### Make Plots ###
    plot.args <- list(
        var = var, object = object, main = NULL, min = min,
        max = max, legend.show = show.legend.plots, data.type = data.type,
        color.panel = color.panel, ...)
    if (!is.null(plot.args$cells.use)) {
        stop("Subsetting with cells.use is incompatible with this function.")
    }
    plots <- .make_vary_cell_plots(
        var, object, vary.cells.meta, cells.meta, vary.cells.levels,
        plot.args, numeric, show.titles)
    # Generate allcells.plot and legend
    plot.args$legend.show <- TRUE
    plot.args$main <- allcells.main
    allcells.plot <- do.call(dittoDimPlot, plot.args)
    legend <- .grab_legend(allcells.plot)
    if (!show.legend.allcells.plot) {
        allcells.plot <- .remove_legend(allcells.plot)
    }
    # Add allcells.plot and legend if wanted
    if (show.allcells.plot) {
        plots$allcells <- allcells.plot
    }
    if (show.legend.single){
        plots$legend <- legend
    }

    #Done!  Return them.
    if(OUT.List){
        return(plots)
    } else {
        if (is.null(ncol) && is.null(nrow)) {
            ncol <- 3
        }
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

.make_vary_cell_plots <- function(
    var, object, vary.cells.meta, cells.meta, vary.cells.levels,
    plot.args, numeric, show.titles) {
    if (numeric) {
        # Case1: Numeric data, colors already matched across plots in min/max
        plots <- lapply(
            vary.cells.levels,
            function(level) {
                if (show.titles) {
                    plot.args$main <- level
                }
                plot.args$cells.use <- cells.meta==level
                do.call(dittoDimPlot, plot.args)
            })
    } else {
        # Case2: non-numeric, colors must be matched accross plots
        all.var.levels <- meta.levels(var, object)
        plots <- lapply(
            vary.cells.levels,
            function(level) {
                if (show.titles) {
                    plot.args$main <- level
                }
                plot.args$cells.use <- cells.meta==level
                plot.levels <-
                    meta.levels(var, object,
                        cells.use = cells.meta==level)
                plot.args$colors <-
                    seq_along(all.var.levels)[all.var.levels %in% plot.levels]
                do.call(dittoDimPlot, plot.args)
            })
    }
    names(plots) <- vary.cells.levels
    plots
}
