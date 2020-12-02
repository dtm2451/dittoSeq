#' Generates dittoDimPlots for multiple features.
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param vars c("var1","var2","var3",...). A vector of vars ('var' in regular \code{\link{dittoDimPlot}}) from which to generate the separate plots.
#' @param ncol,nrow Integer or NULL. How many columns or rows the plots should be arranged into.
#' @param axes.labels.show Logical. Whether axis labels should be shown.
#' Subordinate to \code{xlab} and \code{ylab}.
#' @param list.out Logical. (Default = FALSE) When set to \code{TRUE}, a list of the individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' @param OUT.List Deprecated. Use \code{list.out}
#' @param ...,xlab,ylab,data.out,do.hover,legend.show other parameters passed to \code{\link{dittoDimPlot}}.
#' @return A set of dittoDimPlots either arranged into a grid (default), or output as a list.
#' 
#' @details
#' Given multiple 'var' parameters to \code{vars}, this function creates a \code{\link{dittoDimPlot}} for each one, with minor defaulting tweaks (see below).
#' 
#' By default, these dittoDimPlots are arranged into a grid.
#' Alternatively, if \code{list.out} is set to \code{TRUE}, they are output as a list with each plot named as the \code{vars} being shown.
#' 
#' All parameters that can be adjusted in dittoDimPlot can be adjusted here, but the only input that will change between plots is \code{var}.
#' 
#' @section Slight tweaks to dittoDimPlot defaults:
#' \itemize{
#' \item axes labels are not shown by default to save space (control with \code{axes.labels.show} or \code{xlab} and \code{ylab})
#' \item legends are also not shown to save space (control with \code{legend.show})
#' }
#' 
#' @seealso
#' \code{\link{multi_dittoDimPlotVaryCells}} for an alternate \code{\link{dittoDimPlot}} multi-plotter where the cells/samples are varied between plots.
#' 
#' \code{\link{dittoDimPlot}} for the base dittoDimPlot plotting function and details on all accepted inputs.
#' 
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' 
#' multi_dittoDimPlot(myRNA, c("gene1", "gene2", "clustering"))
#' 
#' # Control grid shape with ncol / nrow
#' multi_dittoDimPlot(myRNA, c("gene1", "gene2", "clustering"),
#'     nrow = 1)
#'     
#' # Output as list instead
#' multi_dittoDimPlot(myRNA, c("gene1", "gene2", "clustering"),
#'     list.out = TRUE)
#'
#' @author Daniel Bunis
#' @export

multi_dittoDimPlot <- function(
    object,
    vars,
    ncol = NULL,
    nrow = NULL,
    axes.labels.show = FALSE,
    list.out = FALSE,
    OUT.List = NULL,
    ...,
    xlab = NA,
    ylab = NA,
    data.out = FALSE,
    do.hover = FALSE,
    legend.show = FALSE
    ) {

    #Interpret axes.labels.show:
    # If xlab/ylab not given, set via axes.labels.show
    lab <- switch(axes.labels.show, "TRUE" = "make", "FALSE" = NULL)
    if (is.na(ylab)) {ylab <- lab}
    if (is.na(xlab)) {xlab <- lab}
    
    if (!is.null(OUT.List)) {
        .Deprecated(msg="Using 'OUT.List', but argument 'OUT.List' is deprecated. Please use 'list.out' instead.")
        list.out <- OUT.List
    }

    if (!list.out && data.out | do.hover) {
        list.out <- TRUE
        message("'data.out' or 'do.hover' requested, outputting as a list.")
    }
    
    if (!length(vars)>=1) {
        stop("No 'vars' provided.")
    }
    
    plots <- lapply(vars, function(X) {
        dittoDimPlot(
            object, X, xlab = xlab, ylab = ylab, data.out = data.out,
            do.hover = do.hover, legend.show = legend.show, ...)
    })
            
    if (list.out){
        names(plots) <- vars
        return(plots)
    } else {
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

#' Generates dittoPlots for multiple features.
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param vars c("var1","var2","var3",...). A vector of gene or metadata names from which to generate the separate plots
#' @param group.by String representing the name of a metadata to use for separating the cells/samples into discrete groups.
#' @param ncol,nrow Integer or NULL. How many columns or rows the plots should be arranged into.
#' @param main,ylab String which sets whether / how plot titles or y-axis labels should be added to each individual plot
#' \itemize{
#' \item When set to \code{"var"}, the \code{vars} names alone will be used.
#' \item When set to \code{"make"}, the default dittoPlot behavior will be observed: For y-axis labels, gene vars will become "'var' expression". Equivalent to "var" for \code{main}.
#' \item When set as any other string, that string will be used as the title / y-axis label for every plot.
#' \item When set to \code{NULL}, titles / axes labels will not be added.
#' }
#' @param OUT.List Deprecated. Use \code{list.out}
#' @param list.out Logical. (Default = FALSE) When set to \code{TRUE}, a list of the individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' @param ...,xlab,data.out,do.hover,legend.show other paramters passed along to \code{\link{dittoPlot}}.
#' @return A set of dittoPlots either arranged into a grid (default), or output as a list.
#' 
#' @details
#' Given multiple 'var' parameters to \code{vars}, this function creates a \code{\link{dittoPlot}} for each one, with minor defaulting tweaks (see below).
#' 
#' By default, these dittoPlots are arranged into a grid.
#' Alternatively, if \code{list.out} is set to \code{TRUE}, they are output as a list with each plot named as the \code{vars} being shown.
#' 
#' All parameters that can be adjusted in dittoPlot can be adjusted here, but the only input that will change between plots is the \code{var}.
#' 
#' @section Slight tweaks to dittoPlot defaults:
#' \itemize{
#' \item axes labels are not shown by default to save space (control with \code{xlab} and \code{ylab})
#' \item legends are also not shown to save space (control with \code{legend.show})
#' }
#' 
#' @seealso
#' \code{\link{dittoPlot}} for the single plot version of this function and details on all accepted inputs. 
#' 
#' \code{\link{dittoDotPlot}} and \code{\link{dittoPlotVarsAcrossGroups}} to show, in a single plot, per-group summaries of the values for multiple vars.
#' 
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' genes <- getGenes(myRNA)[1:4]
#' 
#' multi_dittoPlot(myRNA,
#'     vars = c("gene1", "gene2", "gene3", "gene4"),
#'     group.by = "clustering")
#'
#' #To make it output a grid that is 2x2, to add y-axis labels
#' # instead of titles, and to show legends...
#' multi_dittoPlot(myRNA, c("gene1", "gene2", "gene3", "gene4"), "clustering",
#'     nrow = 2, ncol = 2,           #Make grid 2x2 (only one of these needed)
#'     main = NULL, ylab = "make",   #Add y axis labels instead of titles
#'     legend.show = TRUE)           #Show legends
#'
#' # Output as list instead
#' multi_dittoPlot(myRNA, c("gene1", "gene2", "gene3", "gene4"), "clustering",
#'     list.out = TRUE)
#'
#' @author Daniel Bunis
#' @export

multi_dittoPlot <- function(
    object,
    vars,
    group.by,
    ncol = 3,
    nrow = NULL,
    main = "var",
    ylab = NULL,
    list.out = FALSE,
    OUT.List = NULL,
    ...,
    xlab = NULL,
    data.out = FALSE,
    do.hover = FALSE,
    legend.show = FALSE) {

    if (!is.null(OUT.List)) {
        .Deprecated(msg="Using 'OUT.List', but argument 'OUT.List' is deprecated. Please use 'list.out' instead.")
        list.out <- OUT.List
    }
    
    if (!list.out && data.out | do.hover) {
        list.out <- TRUE
        message("'data.out' or 'do.hover' requested, outputting as a list.")
    }
    
    if (!length(vars)>=1) {
        stop("No 'vars' provided.")
    }
    
    # Prep dittoPlot args
    args <- list(
        object = object, group.by = group.by, xlab = xlab,
        ylab = ylab, main = main, data.out = data.out, do.hover = do.hover,
        legend.show = legend.show, ...)
    
    # Make plots
    plots <- lapply(vars, function(X) {
        
        if (!is.null(ylab)) {
            args$ylab <- ifelse(ylab == "var", X, ylab)
        }
        if (!is.null(main)) {
            args$main <- ifelse(main == "var", X, main)
        }

        do.call(dittoPlot, c(args, var = X))
    })

    #Output
    if (list.out){
        names(plots) <- vars
        return(plots)
    } else {
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

#' Generates multiple dittoDimPlots, for a single feature, where each showing different cells
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed, using colors, for each cell/sample.
#'
#' Alternatively, can be a vector of same length as there are cells/samples in the \code{object}.
#' Discrete or continuous data both work.
#' @param vary.cells.meta String name of a metadata that should be used for selecting which cells to show in each "VaryCells" \code{\link{dittoDimPlot}}.
#' @param vary.cells.levels The values/groupings of the \code{vary.cells.meta} metadata for which to generate a plot.
#' @param show.allcells.plot Logical which sets whether an additional plot showing all of the cells should be added.
#' @param show.legend.plots Logical which sets whether or not legends should be plotted in inidividual VaryCell plots. Default = FALSE.
#' @param show.legend.allcells.plot Logical which sets whether or a legend should be plotted in the allcells plot. Default = FALSE.
#' @param show.legend.single Logical which sets whether to add a single legend as an additional plot. Default = TRUE.
#' @param show.titles Logical which sets whether grouping-levels should be used as titles for the individual VaryCell plots. Default = TRUE.
#' @param ncol,nrow Integers which set dimensions of the plot grid when \code{list.out = FALSE}.
#' @param allcells.main String which adjusts the title of the allcells plot. Default = "All Cells".  Set to \code{NULL} to remove.
#' @param ...,color.panel,colors,min,max,assay,slot,adjustment,data.out,do.hover,swap.rownames additional parameters passed to \code{\link{dittoDimPlot}}.
#' 
#' All parameters of \code{\link{dittoDimPlot}} can be utilized and adjusted except for \code{cells.use}, \code{main}, and \code{legend.show} which are handled with alternative methods here.
#' A few suggestions: \code{reduction.use} for setting which dimensionality reduction space to use.
#' \code{xlab} and \code{ylab} can be set to \code{NULL} to remove the axes labels and provide extra room for the data.
#' \code{size} can be used to adjust the size of the dots.
#' @param OUT.List Deprecated. Use \code{list.out}
#' @param list.out Logical which controls whether the list of plots should be returned as a list instead of as a single grid arrangement of the plots.
#' @return A set of dittoDimPlots either arranged into a grid (default), or output as a list.
#' @details This function generates separate dittoDimPlots that show the same target data, but each for distinct cells.
#' 
#' How cells are separated into distinct plots is controlled with the \code{vary.cells.meta} parameter.
#' Individual \code{\link{dittoDimPlot}}s are created for all levels of \code{var.cells.meta} groupings given to the \code{vary.cells.levels} input (default = all).
#'
#' The function then appends a plot containing all cell/samples when \code{show.allcells.plot = TRUE}, with title of this plot controlled by \code{allcells.main},
#' as well as as single legend when \code{show.legend.single = TRUE}.
#'
#' By default, these dittoDimPlots are output in a grid (default) with \code{ncol} columns and \code{nrow} rows,
#' Alternatively, if \code{list.out} is set to \code{TRUE}, they are returned as a list.
#' In the list, the VaryCell plots will be named by the levels of \code{vary.cells.meta} that they contain,
#' and the optional allcells plot and single legend will be named "allcells" and "legend", respectively.
#'
#' Either continuous or discrete \code{var} data can be displayed.
#' \itemize{
#' \item For continuous data, the range of potential values is calculated at the start, and set, so that colors represent the same value across all plots.
#' \item For discrete data, colors used in each plot are adjusted so that colors represent the same groupings across all plots.
#' }
#'
#' @seealso
#' \code{\link{multi_dittoDimPlot}} for an alternate \code{\link{dittoDimPlot}} multi-plotter where \code{var}s are varied across plots rather than cells/samples
#'
#' \code{\link{dittoDimPlot}} for the base dittoDimPlot plotting function and details on all accepted inputs.
#' 
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' 
#' # This function can be used to quickly scan for differences in expression
#' #   within or across clusters/cell types.
#' multi_dittoDimPlotVaryCells(myRNA, "gene1", vary.cells.meta = "clustering")
#' 
#' # Output as list instead
#' multi_dittoDimPlotVaryCells(myRNA, "gene1", vary.cells.meta = "clustering",
#'     list.out = TRUE)
#'
#' # This function is also great for generating separate plots of each individual
#' #   grouping of a tsne/PCA/umap. This can be useful to check for dispersion
#' #   of groups that might otherwise be hidden behind other cells/samples.
#' #   The effect is similar to faceting, but: all distinct plots are treated
#' #   separately rather than being just a part of the whole, and with portrayal
#' #   of all cells/samples in an additional plot by default.
#' #
#' #   To do so, set 'var' and 'vary.cells.meta' the same.
#' multi_dittoDimPlotVaryCells(myRNA, "clustering", vary.cells.meta = "clustering")
#'
#' # The function can also be used to quickly visualize how separate clustering
#' #   resolutions match up to each other, or perhaps how certain conditions of
#' #   cells disperse across clusters.
#' # (For an alternative method of viewing, and easily quantifying, how discrete
#' #   conditions of cells disperse across clusters, see '?dittoBarPlot')
#' multi_dittoDimPlotVaryCells(myRNA, "groups", vary.cells.meta = "clustering")
#'
#' @author Daniel Bunis
#' @importFrom gridExtra grid.arrange
#' @export

multi_dittoDimPlotVaryCells <- function(
    object,
    var,
    vary.cells.meta,
    vary.cells.levels = metaLevels(vary.cells.meta, object),
    show.titles = TRUE,
    show.allcells.plot = TRUE,
    allcells.main = "All Cells",
    show.legend.single = TRUE,
    show.legend.plots = FALSE,
    show.legend.allcells.plot = FALSE,
    nrow = NULL,
    ncol = NULL,
    list.out = FALSE,
    OUT.List = NULL,
    ...,
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    min = NULL,
    max = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    data.out = FALSE,
    do.hover = FALSE,
    swap.rownames = NULL
    ) {
    
    color.panel <- color.panel[colors]
    
    var.data <- .var_OR_get_meta_or_gene(var, object, assay, slot, adjustment, swap.rownames)
    numeric <- is.numeric(var.data)
    
    vary.data <- meta(vary.cells.meta, object)
    
    if (!is.null(OUT.List)) {
        .Deprecated(msg="Using 'OUT.List', but argument 'OUT.List' is deprecated. Please use 'list.out' instead.")
        list.out <- OUT.List
    }
    
    if (!list.out && data.out | do.hover) {
        list.out <- TRUE
        message("'data.out' or 'do.hover' requested, outputting as a list.")
    }
    
    # Prep plotting args for VaryCell plots
    plot.args <- list(
        object = object, var = var,
        legend.show = show.legend.plots,
        color.panel = color.panel,
        ...,
        assay = assay, slot = slot, adjustment = adjustment,
        min = min, max = max,
        data.out = data.out, do.hover = do.hover,
        swap.rownames = swap.rownames)
    
    if (!is.null(plot.args$cells.use)) {
        stop("Further subsetting with 'cells.use' is not supported.")
    }
    
    if (!is.null(plot.args$main)) {
        message("Universal title adjustment through 'main' ignored. Use 'sub' instead.")
    }
    
    # Make VaryCell plots
    if (numeric) {
        
        plots <- .make_vary_cell_plots_numeric(
            vary.data, vary.cells.levels, var.data,
            min, max,
            plot.args, show.titles)
        
    } else {
        
        plots <- .make_vary_cell_plots_discrete(
            vary.data, vary.cells.levels, var.data,
            all.var.levels = metaLevels(var, object, used.only = FALSE),
            plot.args, show.titles)
    }
    names(plots) <- vary.cells.levels
    
    # Generate and add allcells.plot and legend
    if (show.allcells.plot || show.legend.single) {
        
        plot.args$legend.show <- TRUE
        plot.args$main <- allcells.main
        allcells.plot <- do.call(dittoDimPlot, plot.args)
        legend <- .grab_legend(allcells.plot)
    
        if (show.allcells.plot) {
            
            plots$allcells <- if (show.legend.allcells.plot) {
                allcells.plot
            } else {
                .remove_legend(allcells.plot)
            }
        }
        
        if (show.legend.single){
            plots$legend <- legend
        }
        
    }

    # Done! Return.
    if(list.out) {
        return(plots)
    } else {
        if (is.null(ncol) && is.null(nrow)) {
            ncol <- 3
        }
        return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
    }
}

.make_vary_cell_plots_numeric <- function(
    vary.data, vary.cells.levels,
    var.data, min, max,
    plot.args, show.titles) {
    
    # Set consistent range
    the.range <- range(var.data)
    plot.args$min <- ifelse(is.null(min), the.range[1], min)
    plot.args$max <- ifelse(is.null(max), the.range[2], max)
    
    plots <- lapply(
        vary.cells.levels,
        function(level) {
            
            if (show.titles) {
                plot.args$main <- level
            }
            plot.args$cells.use <- vary.data==level
            
            do.call(dittoDimPlot, plot.args)
        })
}

.make_vary_cell_plots_discrete <- function(
    vary.data, vary.cells.levels,
    var.data, all.var.levels,
    plot.args, show.titles) {
    
    plots <- lapply(
        vary.cells.levels,
        function(level) {
            
            if (show.titles) {
                plot.args$main <- level
            }
            plot.args$cells.use <- vary.data==level
            
            # Ensure consistent colors given that levels may be missing
            plot.levels <- levels(as.factor(var.data[vary.data==level]))
            plot.args$colors <-
                seq_along(all.var.levels)[all.var.levels %in% plot.levels]
            
            do.call(dittoDimPlot, plot.args)
        })
}
