#### multi_dittoDimPlot
#' Generates multiple dittoDimPlots arranged in a grid.
#'
#' @param object A Seurat or SingleCellExperiment object to work with
#' @param vars c("var1","var2","var3",...). A list of vars from which to generate the separate plots
#' @param ncol,nrow Integer/NULL. How many columns or rows the plots should be arranged into
#' @param axes.labels.show Logical. Whether a axis labels should be shown. Ignored if xlab or ylab are set manually.
#' @param OUT.List Logical. (Default = FALSE) When set to \code{TRUE}, a list of the individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' @param legend.show,xlab,ylab,... other paramters passed to \code{\link{dittoDimPlot}}.
#' @return Given multiple 'var' parameters to \code{vars}, this function will output a dittoDimPlot for each one, arranged into a grid, with some slight tweaks to the defaults.
#' If \code{OUT.list} was set to TRUE, the list of individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' All parameters that can be adjusted in dittoDimPlot can be adjusted here, but the only parameter that can be adjusted between each is the \code{var}.
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' genes <- getGenes(myRNA)[1:5]
#' multi_dittoDimPlot(myRNA, c(genes, "clustering"))
#'
#' @author Daniel Bunis
#' @export

multi_dittoDimPlot <- function(
    object,
    vars,
    legend.show = FALSE,
    ncol = NULL,
    nrow = NULL,
    axes.labels.show = FALSE,
    xlab = NA,
    ylab = NA,
    OUT.List = FALSE,
    ...) {

    #Interpret axes.labels.show:
    # If axes.labels.show left as FALSE, set lab to NULL, else "make".
    # Then pass to xlab and ylab unless these were provided.
    lab <- if(!axes.labels.show) {
        NULL
    } else {
        "make"
    }
    if (is.na(ylab)) {ylab <- lab}
    if (is.na(xlab)) {xlab <- lab}

    plots <- lapply(vars, function(X) {
        dittoDimPlot(
            object, X, xlab = xlab, ylab = ylab, legend.show = legend.show, ...)
    })
    if (OUT.List){
        names(plots) <- vars
        return(plots)
    } else {
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

#### multi_dittoPlot
#' Generates multiple dittoPlots arranged into a grid.
#'
#' @param object the Seurat or SingleCellExperiment object to draw from
#' @param vars c("var1","var2","var3",...). A vector of gene or metadata names from which to generate the separate plots
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param color.by String representing the name of a metadata to use for setting color. Default = \code{group.by}.
#' @param ncol,nrow Integers which set how many plots will be arranged per column or per row.
#' Default = 3 columns aand however many rows are required.
#'
#' Set both to NULL to have the grid.arrange function figure out what might be most "square" on its own.
#' @param main,ylab String which sets whether / how plot titles or y-axis labels should be added to each individual plot
#' \itemize{
#' \item When set to \code{"var"}, the \code{vars} names alone will be used.
#' \item When set to \code{"make"}, the default dittoPlot behavior will be observed: Equivalent to "make" for \code{main}, but for y-axis labels, gene vars will become "'var' expression".
#' \item When set as any other string, that string will be used as the title / y-axis label for every plot.
#' \item When set to \code{NULL}, titles / axes labels will not be added.
#' }
#' @param OUT.List Logical. (Default = FALSE) When set to \code{TRUE}, a list of the individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' @param xlab,legend.show,... other paramters passed along to \code{\link{dittoPlot}}.
#' @return Given multiple 'var' parameters, this function will output a dittoPlot for each one, arranged into a grid, just with some slight tweaks to the defaults.
#' If \code{OUT.list} was set to TRUE, the list of individual plots is output instead of the combined multi-plot.
#' All parameters that can be adjusted in dittoPlot can be adjusted here.
#' @seealso
#' \code{\link{dittoPlot}} for the single plot version of this function
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' genes <- getGenes(myRNA)[1:4]
#' multi_dittoPlot(myRNA, genes, group.by = "clustering")
#'
#' # violin-plots in front is often better for large single-cell datasets,
#' # but we cn change the order with 'plots'
#' multi_dittoPlot(myRNA, genes, "clustering",
#'     plots = c("vlnplot","boxplot","jitter"))
#'
#' #To make it output a grid that is 2x2, to add y-axis labels
#' # instead of titles, and to show legends...
#' multi_dittoPlot(myRNA, genes, "clustering",
#'     nrow = 2, ncol = 2,           #Make grid 2x2 (only one of these needed)
#'     main = NULL, ylab = "make",   #Add y axis labels instead of titles
#'     legend.show = TRUE)           #Show legends
#'
#' # We can also facet with 'split.by'
#' multi_dittoPlot(myRNA, genes, "clustering",
#'     split.by = "SNP")
#'
#' @author Daniel Bunis
#' @importFrom ggridges geom_density_ridges2
#' @export

multi_dittoPlot <- function(
    object,
    vars,
    group.by,
    color.by = group.by,
    legend.show = FALSE,
    ncol = 3,
    nrow = NULL,
    main="var",
    ylab = NULL,
    xlab = NULL,
    OUT.List = FALSE,
    ...) {

    plots <- lapply(vars, function(X) {
        args <- list(object, X, group.by, color.by, xlab = xlab,
            ylab = ylab, main = main, legend.show = legend.show, ...)
        if (!is.null(ylab)) {
            args$ylab <- ifelse(ylab == "var", X, ylab)
        }
        if (!is.null(main)) {
            args$main <- ifelse(main == "var", X, main)
        }
        do.call(dittoPlot, args)
    })

    #Output
    if (OUT.List){
        names(plots) <- vars
        return(plots)
    } else {
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}


##################### multi_dittoDimPlotVaryCells #######################
#' Generates multiple dittoDimPlots, each showing different cells, arranged into a grid.
#'
#' @param object A Seurat or SingleCellExperiment object to work with
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample.
#'
#' Alternatively, can be a vector of same length as there are cells/samples in the \code{object}.
#' Discrete or continuous data both work. REQUIRED.
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
#' @param color.panel,colors,min,max,assay,slot,adjustment,... additional parameters passed to \code{\link{dittoDimPlot}}.
#' All parameters except for \code{cells.use}, \code{main}, and \code{legend.show} can be used.
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
#' In the list, the varycells plots will be named by the value of \code{vary.cells.meta} that they contain,
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
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' multi_dittoDimPlotVaryCells(myRNA, "gene1", vary.cells.meta = "clustering")
#'
#' # This function can be used to quickly scan for differences in expression
#' #   within or accross clusters/cell types by providing a gene to 'var'
#' multi_dittoDimPlotVaryCells(myRNA, "gene1", vary.cells.meta = "clustering")
#'
#' # This function is also great for generating separate plots of each individual
#' #   element of a tsne/PCplot/umap. This can be useful to check for dispersion
#' #   of groups that might otherwise be hidden behind other cells/samples.
#' #   To do so, set 'var' and 'vary.cells.meta' the same.
#' multi_dittoDimPlotVaryCells(myRNA, "clustering", vary.cells.meta = "clustering")
#'
#' # The function can also be used to quickly visualize how separate clustering
#' #   resolutions match up to each other, or perhaps how certain conditions of
#' #   cells disperse accross clusters.
#' multi_dittoDimPlotVaryCells(myRNA, "groups", vary.cells.meta = "clustering")
#'
#'
#' # For an alternative method of viewing, and easily quantifying, how discrete
#' #   conditions of cells disperse accross clusters, see '?dittoBarPlot'
#'
#'
#' # Note, for displaying expression or scoring of distinct genes or metadata,
#' #   use 'multi_dittoDimPlot'.  Its split.by variable can then be used to add
#' #   a varyCells-like effect.
#'
#' @author Daniel Bunis
#' @importFrom gridExtra grid.arrange
#' @export

multi_dittoDimPlotVaryCells <- function(
    object, var, vary.cells.meta,
    vary.cells.levels = metaLevels(vary.cells.meta, object),
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL, min = NULL, max = NULL,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    show.titles=TRUE, show.allcells.plot = TRUE, allcells.main = "All Cells",
    show.legend.single = TRUE, show.legend.plots = FALSE,
    show.legend.allcells.plot = FALSE, nrow = NULL, ncol = NULL,
    OUT.List = FALSE, ...)
{
    color.panel <- color.panel[colors]
    cells.meta <- meta(vary.cells.meta,object)
    #Determine if var is continuous vs discrete
    numeric <- ifelse(
        is.numeric(.var_OR_get_meta_or_gene(var, object, assay, slot, adjustment)),
        TRUE,
        FALSE
    )
    if (numeric) {
        the.range <- range(.var_OR_get_meta_or_gene(var, object, assay, slot, adjustment))
        min <- ifelse(is.null(min), the.range[1], min)
        max <- ifelse(is.null(max), the.range[2], max)
    }

    ### Make Plots ###
    plot.args <- list(
        object = object, var = var, main = NULL, min = min,
        max = max, legend.show = show.legend.plots, assay = assay, slot = slot,
        adjustment = adjustment, color.panel = color.panel, ...)
    if (!is.null(plot.args$cells.use)) {
        stop("Further subsetting with 'cells.use' is incompatible with this function.")
    }
    if (!is.null(plot.args$main)) {
        message("Universal title adjustment through 'main' ignored. Use 'sub' instead.")
    }
    plots <- .make_vary_cell_plots(
        object, var, vary.cells.meta, cells.meta, vary.cells.levels,
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
    object, var, vary.cells.meta, cells.meta, vary.cells.levels,
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
        all.var.levels <- metaLevels(var, object)
        plots <- lapply(
            vary.cells.levels,
            function(level) {
                if (show.titles) {
                    plot.args$main <- level
                }
                plot.args$cells.use <- cells.meta==level
                plot.levels <-
                    metaLevels(var, object,
                        cells.use = cells.meta==level)
                plot.args$colors <-
                    seq_along(all.var.levels)[all.var.levels %in% plot.levels]
                do.call(dittoDimPlot, plot.args)
            })
    }
    names(plots) <- vary.cells.levels
    plots
}
