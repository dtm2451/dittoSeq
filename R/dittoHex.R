#' Show RNAseq data, grouped into hexagonal bins, on a scatter or dimensionality reduction plot
#' @name dittoHex
#' 
#' @param bins Numeric or numeric vector giving the number of haxagonal bins in the x and y directions. Set to 30 by default.
#' @param color.method Works differently depending on whether the color.var is continous versus discrete:
#' 
#' \strong{Continuous}: String signifying a function for how target data should be summarized for each bin.
#' Can be any function that summarizes a numeric vector input with a single numeric output value.
#' Default is \code{median}. Other useful options are \code{sum}, \code{mean}, \code{sd}, or \code{mad}.
#' 
#' \strong{Discrete}: A string signifying whether the color should (default) be simply based on the "max" grouping of the bin,
#' or based on the "max.prop"ortion of cells/samples belonging to any grouping.
#' 
#' @param assay,slot,adjustment,assay.x,assay.y,assay.color,assay.extra,slot.x,slot.y,slot.color,slot.extra,adjustment.x,adjustment.y,adjustment.color,adjustment.extra
#' assay, slot, and adjustment set which data to use when the axes, coloring, or \code{extra.vars} are based on expression/counts data. See \code{\link{gene}} for additional information.
#' @param legend.density.title,legend.color.title Strings which set the title for the legends.
#' @param legend.density.breaks,legend.color.breaks Numeric vector which sets the discrete values to label in the density and color.var legends.
#' @param legend.density.breaks.labels,legend.color.breaks.labels String vector, with same length as \code{legend.*.breaks}, which sets the labels for the tick marks of associated legend.
#' @param min.opacity,max.opacity Scalar between [0,1] which sets the minimum or maximum opacities used for the density legend when color is used for \code{color.var} and density is shown via opacity.
#' @param max.color color for highest values of var/max.  Default = blue
#' @param main String, sets the plot title.
#' A default title is automatically generated if based on \code{color.var} and \code{shape.by} when either are provided.
#' To remove, set to \code{NULL}.
#' Defaults to "solid", but see \code{\link[ggplot2]{linetype}} for other options.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot ("p"),
#' a data.frame containing the underlying data for target cells ("data"),
#' and a data.frame containing the underlying data for non-target cells ("Others_data").
#' @inheritParams dittoScatterPlot
#' @inheritParams dittoDimPlot
#' 
#' @return Filler
#' @details
#' This function creates a dataframe with X, Y, color, and faceting data determined by \code{x.var}, \code{y.var}, \code{color.var}, and \code{split.by}.
#' Any extra gene or metadata requested with \code{extra.var} is added as well.
#' For expression/counts data, \code{assay}, \code{slot}, and \code{adjustment} inputs (\code{.x}, \code{.y}, and \code{.color}) can be used to change which data is used, and if it should be adjusted in some way.
#'
#' Next, if a set of cells or samples to use is indicated with the \code{cells.use} input, then the dataframe is split into \code{data} and \code{Others_data} based on subsetting by the target cells/samples.
#'
#' Finally, a hex plot is created using these dataframes.
#' 
#' Fillr filler filler for the rest
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' FILLER
#' \itemize{
#' \item all filler
#' \item \code{size} and \code{opacity} can be used to adjust the size and transparency of the data points.
#' \item Colors used can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes used can be adjusted with \code{shape.panel}.
#' \item Color and shape labels can be changed using \code{rename.color.groups} and \code{rename.shape.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab completion lookup.
#' }
#'
#' @seealso
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{x.var}, \code{y.var}, \code{color.var}, \code{shape.by}, and \code{hover.data} options are.
#'
#' \code{\link{dittoDimPlot}} for making very similar data representations, but where dimensionality reduction (PCA, t-SNE, UMAP, etc.) dimensions are the scatterplot axes.
#'
#' @author Daniel Bunis with some code adapted from Giuseppe D'Agostino
#' @examples
#' # dittoSeq handles bulk and single-cell data quit similarly.
#' # The SingleCellExperiment object structure is used for both,
#' # but all functions can be used similarly directly on Seurat
#' # objects as well.
#'
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#' 
#' # Mock up some nCount_RNA and nFeature_RNA metadata
#' #  == the default way to extract
#' myRNA$nCount_RNA <- runif(60,200,1000)
#' myRNA$nFeature_RNA <- myRNA$nCount_RNA*runif(60,0.95,1.05)
#' # and also percent.mito metadata
#' myRNA$percent.mito <- sample(c(runif(50,0,0.05),runif(10,0.05,0.2)))
#'
#' dittoScatterHex(
#'     myRNA, x.var = "nCount_RNA", y.var = "nFeature_RNA")
#' dittoDimHex(myRNA)
#' 
#' # We don't have too many samples here, so let's increase the bin size.
#' dittoDimHex(myRNA, bins = 10)
#' 
#' # x and y bins can be set separately, useful for non-square plots
#' dittoDimHex(myRNA, bins = c(20, 10))
#' 
#' ### Coloring
#' # Default coloring, as above, is by cell/sample density in the region, but
#' # 'color.var' can be used to color the data by another metric.
#' # Density with then be represented via bin opacity.
#' dittoDimHex(myRNA, color.var = "groups", bins = 10)
#' dittoDimHex(myRNA, color.var = "gene1", bins = 10)
#' 
#' # 'color.method' is then used to adjust how the target data is summarized
#' dittoDimHex(myRNA, color.var = "groups", bins = 10,
#'     color.method = "max.prop")
#' dittoDimHex(myRNA, color.var = "gene1", bins = 10,
#'     color.method = "mean")
#' 
#' ### Additional Features:
#' 
#' # Faceting with 'split.by'
#' dittoDimHex(myRNA, bins = 10, split.by = "groups")
#' dittoDimHex(myRNA, bins = 10, split.by = c("groups", "clustering"))
#' 
#' # Underlying data output with 'data.out = TRUE'
#' dittoDimHex(myRNA, data.out = TRUE)
#' 
#' # Contour lines can be added with 'do.contours = TRUE'
#' dittoDimHex(myRNA, bins = 10,
#'     do.contour = TRUE,
#'     contour.color = "lightblue", # Optional, black by default
#'     contour.linetype = "dashed") # Optional, solid by default
#' 
#' # Trajectories can be added to dittoDimHex plots (see above for details)
#' dittoDimHex(myRNA, bins = 10,
#'     add.trajectory.lineages = list(c(1,2,4), c(1,4), c(1,3)),
#'     trajectory.cluster.meta = "clustering")
NULL

#' @describeIn dittoHex Make a scatter plot of RNAseq data, grouped into hexagonal bins
#' @export
dittoScatterHex <- function(
    object, x.var, y.var, color.var = NULL, bins = 30,
    color.method = "median",
    split.by = NULL, extra.vars = NULL,
    cells.use = NULL,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    split.nrow = NULL,
    split.ncol = NULL,
    assay.x = .default_assay(object),
    slot.x = .default_slot(object),
    adjustment.x = NULL,
    assay.y = .default_assay(object),
    slot.y = .default_slot(object),
    adjustment.y = NULL,
    assay.color = .default_assay(object),
    slot.color = .default_slot(object),
    adjustment.color = NULL,
    assay.extra = .default_assay(object),
    slot.extra = .default_slot(object),
    adjustment.extra = NULL,
    min.opacity = 0.2,
    max.opacity = 1,
    min.color = "#F0E442",
    max.color = "#0072B2",
    min = NULL,
    max = NULL,
    xlab = x.var,
    ylab = y.var,
    main = "make",
    sub = NULL,
    theme = theme_bw(),
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    legend.show = TRUE,
    legend.color.title = "make",
    legend.color.size = 5,
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.density.title = if (isBulk(object)) "Samples" else "Cells",
    legend.density.breaks = waiver(),
    legend.density.breaks.labels = waiver(),
    data.out = FALSE) {

    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)

    # Make dataframe
    data <- .scatter_data_gather(
        object = object, cells.use = cells.use, x.var = x.var, y.var = y.var,
        color.var = color.var, shape.by = NULL,
        split.by = split.by, extra.vars = extra.vars,
        assay.x = assay.x, slot.x = slot.x, adjustment.x = adjustment.x,
        assay.y = assay.y, slot.y = slot.y, adjustment.y = adjustment.y,
        assay.color = assay.color, slot.color = slot.color,
        adjustment.color = adjustment.color,
        assay.extra = assay.extra, slot.extra = slot.extra,
        adjustment.extra = adjustment.extra
    )[cells.use,]

    # Set titles if "make"
    main <- .leave_default_or_null(main,
        "filler")
    legend.color.title <- .leave_default_or_null(legend.color.title,
        paste(color.var, color.method, sep = ",\n"))

    # Make the plot
    p <- .ditto_scatter_hex(
        data, bins, color.method, color.panel, colors,
        min.opacity, max.opacity, min.color, max.color, min, max,
        xlab, ylab, main, sub, theme, legend.show, legend.color.title,
        legend.color.size, legend.color.breaks, legend.color.breaks.labels,
        legend.density.title, legend.density.breaks, legend.density.breaks.labels)
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, object, cells.use)
    }
    
    ### Add extra features
    if (do.contour) {
        p <- .add_contours(p, data, contour.color,  contour.linetype)
    }

    ### RETURN the PLOT ###
    if (data.out) {
        return(list(plot = p, data = data))
    } else{
        return(p)
    }
}

#' @describeIn dittoHex Show RNAseq data overlayed on a tsne, pca, or similar, grouped into hexagonal bins
#' @export
dittoDimHex <- function(
    object, color.var = NULL, bins = 30,
    color.method = "median",
    reduction.use = .default_reduction(object), dim.1 = 1, dim.2 = 2,
    cells.use = NULL,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    split.by = NULL, extra.vars = NULL,
    split.nrow = NULL, split.ncol = NULL,
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL,
    assay.extra = assay, slot.extra = slot, adjustment.extra = adjustment,
    show.axes.numbers = TRUE,
    show.grid.lines = !grepl("umap|tsne", tolower(reduction.use)),
    main = "make", sub = NULL, xlab = "make", ylab = "make",
    theme = theme_bw(),
    do.contour = FALSE, contour.color = "black", contour.linetype = 1,
    min.opacity = 0.2, max.opacity = 1,
    min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
    add.trajectory.lineages = NULL, add.trajectory.curves = NULL,
    trajectory.cluster.meta, trajectory.arrow.size = 0.15, data.out = FALSE,
    legend.show = TRUE,
    legend.color.title = "make",
    legend.color.size = 5,
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.density.title = if (isBulk(object)) "Samples" else "Cells",
    legend.density.breaks = waiver(),
    legend.density.breaks.labels = waiver()
    ) {

    # Generate the x/y dimensional reduction data and plot titles.
    xdat <- .extract_Reduced_Dim(reduction.use, dim.1, object)
    ydat <- .extract_Reduced_Dim(reduction.use, dim.2, object)
    xlab <- .leave_default_or_null(xlab, xdat$name)
    ylab <- .leave_default_or_null(ylab, ydat$name)
    main <- .leave_default_or_null(main, color.var, length(color.var)!=1)

    # Edit theme
    if (!show.grid.lines) {
        theme <- theme + theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    }
    if (!show.axes.numbers) {
        theme <- theme +
            theme(axis.text.x=element_blank(), axis.text.y=element_blank())
    }

    # Make dataframes and plot
    p.df <- dittoScatterHex(
        object, xdat$embeddings, ydat$embeddings, color.var, bins,
        color.method, split.by,
        extra.vars, cells.use, color.panel, colors,
        split.nrow, split.ncol, NA, NA, NA, NA, NA, NA,
        assay, slot, adjustment, assay.extra, slot.extra, adjustment.extra,
        min.opacity, max.opacity, min.color, max.color, min, max,
        xlab, ylab, main, sub, theme,
        do.contour, contour.color, contour.linetype,
        legend.show, legend.color.title, legend.color.size,
        legend.color.breaks, legend.color.breaks.labels, legend.density.title,
        legend.density.breaks, legend.density.breaks.labels, data.out = TRUE)
    p <- p.df$plot
    data <- p.df$data

    # Add extra features
    if (is.list(add.trajectory.lineages)) {
        p <- .add_trajectory_lineages(
            p, add.trajectory.lineages, trajectory.cluster.meta,
            trajectory.arrow.size, object, reduction.use, dim.1, dim.2)
    }

    if (is.list(add.trajectory.curves)) {
        p <- .add_trajectory_curves(
            p, add.trajectory.curves, trajectory.arrow.size, dim.1, dim.2)
    }
    
    ### RETURN the PLOT ###
    if (data.out) {
        return(list(
            plot = p,
            data = data))
    } else {
        return(p)
    }
}

.ditto_scatter_hex <- function(
    data,
    bins,
    color.method,
    color.panel,
    colors,
    min.opacity,
    max.opacity,
    min.color,
    max.color,
    min,
    max,
    xlab,
    ylab,
    main,
    sub,
    theme,
    legend.show,
    legend.color.title,
    legend.color.size,
    legend.color.breaks,
    legend.color.breaks.labels,
    legend.density.title,
    legend.density.breaks,
    legend.density.breaks.labels
) {

    ### Set up plotting
    p <- ggplot() + ylab(ylab) + xlab(xlab) + ggtitle(main,sub) + theme

    ### Determine how to add data while adding proper theming
    aes.args <- list(x = "X", y = "Y")
    geom.args <- list(
        data = data, bins = bins, na.rm = TRUE)
    
    if (is.null(data$color)) {
        ## Simple case: no color.var given
        
        # Set color scale based on density
        p <- p + scale_fill_gradient(
            name = legend.density.title,
            low= min.color,
            high = max.color,
            breaks = legend.density.breaks,
            labels = legend.density.breaks.labels)
        
    } else {
        ## Setup for ggplot.multistats::stat_summaries_hex
        .error_if_no_ggplot.multistats()
        
        # Set alpha scale based on density
        p <- p + scale_alpha_continuous(
            name = legend.density.title,
            range = c(min.opacity, max.opacity),
            breaks = legend.density.breaks,
            labels = legend.density.breaks.labels)
        
        # Set aesthetics
        aes.args$z <- "color"
        aes.args$fill <- "stat(c)"
        aes.args$alpha <- "stat(d)"
        # Fix for when color metadata is a factor
        aes.args$group <- 1
        
        # Determine how 'c' and 'd' should be calculated &
        # set fill based on color.var and the color.method
        if (is.numeric(data$color) || is.character(color.method) && color.method == "max.prop") {
            
            geom.args$funs <- c(
                c = if (is.numeric(data$color)) {
                    color.method
                } else {
                    function(x) max(table(x)/length(x))
                }, d = length)
            
            p <- p + scale_fill_gradient(
                name = legend.color.title,
                low= min.color,
                high = max.color,
                limits = c(
                    ifelse(is.null(min), NA, min),
                    ifelse(is.null(max), NA, max)),
                breaks = legend.color.breaks,
                labels = legend.color.breaks.labels)
        
        } else {
            
            geom.args$funs <- 
                c(c = function(x) names(which.max(table(x))), d = length)
            
            # Fix title in default case where color.method was left as 'median'
            legend.color.title <- gsub("\nmedian", "\nmax", legend.color.title)
            
            p <- p + scale_fill_manual(
                    name = legend.color.title,
                    values = color.panel[colors]) +
                # Set legend key size
                guides(color = guide_legend(override.aes = list(size=legend.color.size)))
        }
    }
    
    ### Add data
    geom.args$mapping <- do.call(aes_string, aes.args)
    if (!is.null(data$color)) {
        p <- p + do.call(ggplot.multistats::stat_summaries_hex, geom.args)
    } else {
        p <- p + do.call(stat_bin_hex, geom.args)
    }

    if (!legend.show) {
        p <- .remove_legend(p)
    }

    p
}

.scatter_data_gather <- function(
    object,
    cells.use,
    x.var,
    y.var,
    color.var,
    shape.by,
    split.by,
    extra.vars,
    assay.x,
    slot.x,
    adjustment.x,
    assay.y,
    slot.y,
    adjustment.y,
    assay.color,
    slot.color,
    adjustment.color,
    assay.extra,
    slot.extra,
    adjustment.extra,
    do.hover = FALSE,
    hover.data = NULL,
    hover.assay = NULL,
    hover.slot = NULL,
    hover.adjustment = NULL,
    rename.color.groups = NULL,
    rename.shape.groups = NULL
    ) {

    all.cells <- .all_cells(object)
    
    # Make dataframe
    vars <- list(x.var, y.var, color.var, shape.by)
    names <- list("X", "Y", "color", "shape")
    assays <- list(assay.x, assay.y, assay.color, NA)
    slots <- list(slot.x, slot.y, slot.color, NA)
    adjustments <- list(adjustment.x, adjustment.y, adjustment.color, NA)
    relabels <- list(NULL, NULL, rename.color.groups, rename.shape.groups)

    dat <- data.frame(row.names = all.cells)
    for (i in seq_along(vars)) {
        dat <- .add_by_cell(dat, vars[[i]], names[[i]], object, assays[[i]],
            slots[[i]], adjustments[[i]], NULL, relabels[[i]])
    }

    extra.vars <- unique(c(split.by, extra.vars))
    dat <- .add_by_cell(dat, extra.vars, extra.vars, object, assay.extra,
        slot.extra, adjustment.extra, mult = TRUE)

    if (do.hover) {
        dat$hover.string <- .make_hover_strings_from_vars(
            hover.data, object, hover.assay, hover.slot, hover.adjustment)
    }
    
    dat
}
