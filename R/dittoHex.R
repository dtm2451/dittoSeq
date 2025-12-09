#' Show RNAseq data, grouped into hexagonal bins, on a scatter or dimensionality reduction plot
#' @name dittoHex
#' 
#' @param x.var,y.var Single string giving a gene or metadata that will be used for the x- and y-axis of the scatterplot.
#' Note: must be continuous.
#'
#' Alternatively, can be a directly supplied numeric vector of length equal to the total number of cells/samples in \code{object}.
#' 
#' @param rename.color.groups String vector containing new names for the identities of discrete color groups.
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}.
#' @param xlab,ylab Strings which set the labels for the axes. To remove, set to \code{NULL}.
#' @param bins Numeric or numeric vector giving the number of haxagonal bins in the x and y directions. Set to 30 by default.
#' @param color.method Works differently depending on whether the \code{color.var}-data is continuous versus discrete:
#' 
#' \strong{Continuous}: String naming a function for how target data should be summarized for each bin.
#' Can be any function that summarizes a numeric vector input with a single numeric output value.
#' Default is \code{median}. Other useful options are \code{sum}, \code{mean}, \code{sd}, or \code{mad}.
#' You can also use a previously assigned function; e.g. first run \code{pNonZero <- function(x) \{ sum(x!=0)/length(x)\}},
#' then you give \code{color.method = "pNonZero"}
#' 
#' \strong{Discrete}: A string signifying whether the color should (default) be simply based on the "max" grouping of the bin,
#' based on "prop.<value>" the proportion of a specific value (e.g. "prop.A" or "prop.TRUE"),
#' or based on the "max.prop"ortion of cells/samples belonging to any grouping.
#'
#' @param legend.density.title,legend.color.title Strings which set the title for the legends.
#' @param legend.density.breaks,legend.color.breaks Numeric vector which sets the discrete values to label in the density and color.var legends.
#' @param legend.density.breaks.labels,legend.color.breaks.labels String vector, with same length as \code{legend.*.breaks}, which sets the labels for the tick marks or hex icons of the associated legend.
#' @param min.opacity,max.opacity Scalar between [0,1] which sets the minimum or maximum opacity used for the density legend (when color is used for \code{color.var} data and density is shown via opacity).
#' @param min.density,max.density Number which sets the min/max values used for the density scale.
#' Used no matter whether density is represented through opacity or color.
#' @param min.color,max.color color for the min/max values of the color scale. 
#' @param min,max Number which sets the values associated with the minimum or maximum color for \code{color.var} data.
#' @param main String, sets the plot title. The default title is either "Density", \code{color.var}, or NULL, depending on the identity of \code{color.var}.
#' To remove, set to \code{NULL}.
#' @param data.out Logical. When set to \code{TRUE}, changes the output from the plot alone to a list containing the plot ("plot"),
#' and data.frame of the underlying data for target cells ("data").
#' @param add.trajectory.curves List of matrices, each representing coordinates for a trajectory path, from start to end, where matrix columns represent x (\code{dim.1}) and y (\code{dim.2}) coordinates of the paths.
#'
#' Alternatively, (for dittoDimHex only, but not dittoScatterHex) a list of lists(/princurve objects) can be provided.
#' Thus, if the \code{\link[slingshot]{slingshot}} package was used for trajectory analysis,
#' you can provide \code{add.trajectory.curves = slingCurves('object')}
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a named list containing:\itemize{
#' \item "p": the plot
#' \item "data": a data.frame containing the underlying data
#' \item "cols_used": a named list providing the columns of 'data' ultimately used in plotting the named elements
#' \item "to_dittoViz": the dataframe extracted by dittoSeq and passed to \code{dittoViz::\link[dittoViz]{scatterHex}} for plotting.
#' }
#' @inheritParams dittoScatterPlot
#' @inheritParams dittoDimPlot
#' 
#' @details
#' The functions create a dataframe with x and y coordinates for each cell/sample, determined by either \code{x.var} and \code{y.var} for \code{dittoScatterHex},
#' or \code{reduction.use}, \code{dim.1} (x), and \code{dim.2} (y) for \code{dittoDimHex}.
#' Extra data requested by \code{color.var} for coloring, \code{split.by} for faceting, or \code{extra.var} for manual external manipulations, are added to the dataframe as well.
#' For expression/counts data, \code{assay}, \code{slot}, and \code{adjustment} inputs can be used to select which values to use, and if they should be adjusted in some way.
#'
#' The dataframe is then subset to only target cells/samples based on the \code{cells.use} input.
#'
#' Finally, a hex plot is created using this dataframe:
#' 
#' If \code{color.var} is not rovided, coloring is based on the density of cells/samples within each hex bin.
#' When \code{color.var} is provided, density is represented through opacity while coloring is based on a summarization, chosen with the \code{color.method} input, of the target \code{color.var} data.
#' 
#' If \code{split.by} was used, the plot will be split into a matrix of panels based on the associated groupings.
#'
#' @return A ggplot object where colored hexagonal bins are used to summarize RNAseq data in a scatterplot or tSNE, PCA, UMAP.
#'
#' Alternatively, if \code{data.out=TRUE}, a named list containing four elements. See the description of that argument above for further details.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item Colors: \code{min.color} and \code{max.color} adjust the colors for continuous data.
#' \item For discrete \code{color.var} plotting with \code{color.method = "max"}, colors are instead adjusted with \code{color.panel} and/or \code{colors} & the labels of the groupings can be changed using \code{rename.color.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.color.title} and \code{legend.density.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab completion lookup.
#' }
#' 
#' @section Additional Features:
#' Other tweaks and features can be added as well.
#' Each is accessible through 'tab' autocompletion starting with "\code{do.}"\code{---} or "\code{add.}"\code{---},
#' and if additional inputs are involved in implementing or tweaking these, the associated inputs will start with the "\code{---.}":
#' \itemize{
#' \item If \code{do.contour} is provided, density gradiant contour lines will be overlaid with color and linetype adjustable via \code{contour.color} and \code{contour.linetype}.
#' \item If \code{add.trajectory.lineages} is provided a list of vectors (each vector being cluster names from start-cluster-name to end-cluster-name), and a metadata name pointing to the relevant clustering information is provided to \code{trajectory.cluster.meta},
#' then median centers of the clusters will be calculated and arrows will be overlayed to show trajectory inference paths in the current dimmenionality reduction space.
#' \item If \code{add.trajectory.curves} is provided a list of matrices (each matrix containing x, y coordinates from start to end), paths and arrows will be overlayed to show trajectory inference curves in the current dimmenionality reduction space.
#' Arrow size is controlled with the \code{trajectory.arrow.size} input.
#' }
#'
#' @seealso
#' \code{\link{dittoDimPlot}} and \code{\link{dittoScatterPlot}} for making very similar data representations, but where each cell is represented individually.
#' It is often best to investigate your data with both the individual and hex-bin methods, then pick whichever is the best representation for your particular goal.
#' 
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{var}, \code{split.by}, etc. options are of an \code{object}.
#'
#' \code{\link{getReductions}} to see what the \code{reduction.use} options are of an \code{object}.
#' 
#' @author Daniel Bunis with some code adapted from Giuseppe D'Agostino
#' @examples
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
#' dittoDimHex(myRNA, color.var = "clustering", bins = 10)
#' dittoDimHex(myRNA, color.var = "gene1", bins = 10)
#' 
#' # 'color.method' is then used to adjust how the target data is summarized
#' dittoDimHex(myRNA, color.var = "groups", bins = 10,
#'     color.method = "max.prop")
#' dittoDimHex(myRNA, color.var = "gene1", bins = 10,
#'     color.method = "mean")
#' # One particularly useful 'color.method' for discrete 'color.var'-data is
#' #   to use 'prop.<value>' to color by the proportion of a particular value
#' #   within each bin:
#' dittoDimHex(myRNA, color.var = "groups", bins = 10,
#'     color.method = "prop.A")
#' 
#' ### Additional Features:
#' 
#' # Faceting with 'split.by'
#' dittoDimHex(myRNA, bins = 10, split.by = "groups")
#' dittoDimHex(myRNA, bins = 10, split.by = c("groups", "clustering"))
#' 
#' # Faceting can also be used to show multiple continuous variables side-by-side
#' #   by giving a vector of continuous metadata or gene names to 'color.var'.
#' #   This can also be combined with 1 'split.by' variable, with direction then
#' #   controlled via 'multivar.split.dir':
#' dittoDimHex(myRNA, bins = 10,
#'     color.var = c("gene1", "gene2"))
#' dittoDimHex(myRNA, bins = 10,
#'     color.var = c("gene1", "gene2"),
#'     split.by = "groups")
#' dittoDimHex(myRNA, bins = 10,
#'     color.var = c("gene1", "gene2"),
#'     split.by = "groups",
#'     multivar.split.dir = "row")
#' 
#' # Underlying data output with 'data.out = TRUE'
#' dittoDimHex(myRNA, data.out = TRUE)
#' 
#' # Contour lines can be added with 'do.contours = TRUE'
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'     dittoDimHex(myRNA, bins = 10,
#'         do.contour = TRUE,
#'         contour.color = "lightblue", # Optional, black by default
#'         contour.linetype = "dashed") # Optional, solid by default
#' }
#' 
#' # Trajectories can be added to dittoDimHex plots
#' dittoDimHex(myRNA, bins = 10,
#'     add.trajectory.lineages = list(c(1,2,4), c(1,4), c(1,3)),
#'     trajectory.cluster.meta = "clustering")
NULL

#' @describeIn dittoHex Show RNAseq data overlayed on a tsne, pca, or similar, grouped into hexagonal bins
#' @export
dittoDimHex <- function(
    object,
    color.var = NULL,
    bins = 30,
    color.method = NULL,
    reduction.use = .default_reduction(object),
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    split.by = NULL,
    extra.vars = NULL,
    multivar.split.dir = c("col", "row"),
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL,
    assay.extra = assay,
    slot.extra = slot,
    adjustment.extra = adjustment,
    show.axes.numbers = TRUE,
    show.grid.lines = if (is.character(reduction.use)) { !grepl("umap|tsne", tolower(reduction.use)) } else {TRUE},
    main = "make",
    sub = NULL,
    xlab = "make",
    ylab = "make",
    theme = theme_bw(),
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    min.density = NA,
    max.density = NA,
    min.color = "#F0E442",
    max.color = "#0072B2",
    min.opacity = 0.2,
    max.opacity = 1, 
    min = NA,
    max = NA,
    rename.color.groups = NULL,
    do.ellipse = FALSE,
    do.label = FALSE,
    labels.size = 5,
    labels.highlight = TRUE,
    labels.use.numbers = FALSE,
    labels.numbers.spacer = ": ",
    labels.repel = TRUE,
    labels.split.by = split.by,
    labels.repel.adjust = list(),
    add.trajectory.lineages = NULL,
    add.trajectory.curves = NULL,
    trajectory.cluster.meta = NULL,
    trajectory.arrow.size = 0.15,
    add.xline = NULL,
    xline.linetype = "dashed",
    xline.color = "black",
    xline.linewidth = 0.5,
    xline.opacity = 1,
    add.yline = NULL,
    yline.linetype = "dashed",
    yline.color = "black",
    yline.linewidth = 0.5,
    yline.opacity = 1,
    add.abline = NULL,
    abline.slope = 1,
    abline.linetype = "solid",
    abline.color = "black",
    abline.linewidth = 0.5,
    abline.opacity = 1,
    data.out = FALSE,
    legend.show = TRUE,
    legend.color.title = "make",
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.density.title = if (isBulk(object)) "Samples" else "Cells",
    legend.density.breaks = waiver(),
    legend.density.breaks.labels = waiver()
    ) {
    
    multivar.split.dir <- match.arg(multivar.split.dir)

    # Generate the x/y dimensional reduction data and plot titles.
    xdat <- .extract_Reduced_Dim(reduction.use, dim.1, object)
    ydat <- .extract_Reduced_Dim(reduction.use, dim.2, object)
    xlab <- .leave_default_or_null(xlab, xdat$name)
    ylab <- .leave_default_or_null(ylab, ydat$name)

    # Make dataframes and plot
    dittoScatterHex(
        object,
        x.var = xdat$embeddings,
        y.var = ydat$embeddings,
        color.var = color.var,
        bins = bins,
        color.method = color.method,
        split.by = split.by,
        extra.vars = extra.vars,
        cells.use = cells.use,
        color.panel = color.panel,
        colors = colors,
        multivar.split.dir = multivar.split.dir,
        split.nrow = split.nrow,
        split.ncol = split.ncol,
        split.adjust = split.adjust,
        assay.x = NA,
        slot.x = NA,
        adjustment.x = NULL,
        assay.y = NA,
        slot.y = NA,
        adjustment.y = NULL,
        assay.color = assay,
        slot.color = slot,
        adjustment.color = adjustment,
        assay.extra = assay.extra,
        slot.extra = slot.extra,
        adjustment.extra = adjustment.extra,
        swap.rownames = swap.rownames,
        min.density = min.density,
        max.density = max.density,
        min.color = min.color,
        max.color = max.color,
        min.opacity = min.opacity,
        max.opacity = max.opacity,
        min = min,
        max = max,
        rename.color.groups = rename.color.groups,
        show.grid.lines = show.grid.lines,
        show.axes.numbers = show.axes.numbers,
        xlab = xlab,
        ylab = ylab,
        main = main,
        sub = sub,
        theme = theme,
        do.contour = do.contour,
        contour.color = contour.color,
        contour.linetype = contour.linetype,
        do.ellipse = do.ellipse,
        do.label = do.label,
        labels.size = labels.size,
        labels.highlight = labels.highlight,
        labels.use.numbers = labels.use.numbers,
        labels.numbers.spacer = labels.numbers.spacer,
        labels.repel = labels.repel,
        labels.split.by = labels.split.by,
        labels.repel.adjust = labels.repel.adjust,
        add.trajectory.lineages = add.trajectory.lineages,
        add.trajectory.curves = add.trajectory.curves,
        trajectory.cluster.meta = trajectory.cluster.meta,
        trajectory.arrow.size = trajectory.arrow.size,
        add.xline = add.xline,
        xline.linetype = xline.linetype,
        xline.color = xline.color,
        xline.linewidth = xline.linewidth,
        xline.opacity = xline.opacity,
        add.yline = add.yline,
        yline.linetype = yline.linetype,
        yline.color = yline.color,
        yline.linewidth = yline.linewidth,
        yline.opacity = yline.opacity,
        add.abline = add.abline,
        abline.slope = abline.slope,
        abline.linetype = abline.linetype,
        abline.color = abline.color,
        abline.linewidth = abline.linewidth,
        abline.opacity = abline.opacity,
        legend.show = legend.show,
        legend.color.title = legend.color.title,
        legend.color.breaks = legend.color.breaks,
        legend.color.breaks.labels = legend.color.breaks.labels,
        legend.density.title = legend.density.title,
        legend.density.breaks = legend.density.breaks,
        legend.density.breaks.labels = legend.density.breaks.labels,
        data.out = data.out)
}

#' @describeIn dittoHex Make a scatter plot of RNAseq data, grouped into hexagonal bins
#' @export
dittoScatterHex <- function(
    object,
    x.var,
    y.var,
    color.var = NULL,
    bins = 30,
    color.method = NULL,
    split.by = NULL,
    extra.vars = NULL,
    cells.use = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    multivar.split.dir = c("col", "row"),
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
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
    swap.rownames = NULL,
    min.density = NA,
    max.density = NA,
    min.color = "#F0E442",
    max.color = "#0072B2",
    min.opacity = 0.2,
    max.opacity = 1,
    min = NA,
    max = NA,
    rename.color.groups = NULL,
    show.grid.lines = TRUE,
    show.axes.numbers = TRUE,
    xlab = x.var,
    ylab = y.var,
    main = "make",
    sub = NULL,
    theme = theme_bw(),
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    do.ellipse = FALSE,
    do.label = FALSE,
    labels.size = 5,
    labels.highlight = TRUE,
    labels.use.numbers = FALSE,
    labels.numbers.spacer = ": ",
    labels.repel = TRUE,
    labels.split.by = split.by,
    labels.repel.adjust = list(),
    add.trajectory.lineages = NULL,
    add.trajectory.curves = NULL,
    trajectory.cluster.meta = NULL,
    trajectory.arrow.size = 0.15,
    add.xline = NULL,
    xline.linetype = "dashed",
    xline.color = "black",
    xline.linewidth = 0.5,
    xline.opacity = 1,
    add.yline = NULL,
    yline.linetype = "dashed",
    yline.color = "black",
    yline.linewidth = 0.5,
    yline.opacity = 1,
    add.abline = NULL,
    abline.slope = 1,
    abline.linetype = "solid",
    abline.color = "black",
    abline.linewidth = 0.5,
    abline.opacity = 1,
    legend.show = TRUE,
    legend.color.title = "make",
    legend.color.breaks = waiver(),
    legend.color.breaks.labels = waiver(),
    legend.density.title = if (isBulk(object)) "Samples" else "Cells",
    legend.density.breaks = waiver(),
    legend.density.breaks.labels = waiver(),
    data.out = FALSE) {

    # Standardize cells/samples vectors.
    cells.use <- .which_cells(cells.use, object)
    multivar.split.dir <- match.arg(multivar.split.dir)

    # Make dataframe
    pulled_data <- .data_gather_to_df(
        object, color.var,
        x.by = x.var, y.by = y.var,
        extra.vars = unique(c(split.by, extra.vars, trajectory.cluster.meta)),
        assay = assay.color, slot = slot.color,
        swap.rownames = swap.rownames,
        x.assay = assay.x, x.slot = slot.x,
        y.assay = assay.y, y.slot = slot.y,
        extra.assay = assay.extra, extra.slot = slot.extra
    )
    
    if (!show.axes.numbers) {
        theme <- theme +
            theme(axis.text.x=element_blank(), axis.text.y=element_blank())
    }

    # Make the plot
    viz_out <- dittoViz::scatterHex(
        data_frame = pulled_data,
        x.by = if ('_x.by' %in% colnames(pulled_data)) {'_x.by'} else {x.var},
        y.by = if ('_y.by' %in% colnames(pulled_data)) {'_y.by'} else {y.var},
        color.by = if ('_var' %in% colnames(pulled_data)) {'_var'} else {color.var},
        bins = bins,
        color.method = color.method,
        split.by = split.by,
        rows.use = cells.use,
        color.panel = color.panel,
        colors = colors,
        x.adjustment = adjustment.x,
        y.adjustment = adjustment.y,
        color.adjustment = adjustment.color,
        x.adj.fxn = NULL,
        y.adj.fxn = NULL,
        color.adj.fxn = NULL,
        multivar.split.dir = multivar.split.dir,
        split.nrow = split.nrow,
        split.ncol = split.ncol,
        split.adjust = split.adjust,
        min.density = min.density,
        max.density = max.density,
        min.color = min.color,
        max.color = max.color,
        min.opacity = min.opacity,
        max.opacity = min.opacity,
        min = min,
        max = max,
        rename.color.groups = rename.color.groups,
        xlab = xlab,
        ylab = ylab,
        main = main,
        sub = sub,
        theme = theme,
        do.contour = do.contour,
        contour.color = contour.color,
        contour.linetype = contour.linetype,
        do.ellipse = do.ellipse,
        do.label = do.label,
        labels.size = labels.size,
        labels.highlight = labels.highlight,
        labels.use.numbers = labels.use.numbers,
        labels.numbers.spacer = labels.numbers.spacer,
        labels.repel = labels.repel,
        labels.split.by = labels.split.by,
        labels.repel.adjust = labels.repel.adjust,
        add.trajectory.by.groups = add.trajectory.lineages,
        add.trajectory.curves = add.trajectory.curves,
        trajectory.group.by = trajectory.cluster.meta,
        trajectory.arrow.size = trajectory.arrow.size,
        add.xline = add.xline,
        xline.linetype = xline.linetype,
        xline.color = xline.color,
        xline.linewidth = xline.linewidth,
        xline.opacity = xline.opacity,
        add.yline = add.yline,
        yline.linetype = yline.linetype,
        yline.color = yline.color,
        yline.linewidth = yline.linewidth,
        yline.opacity = yline.opacity,
        add.abline = add.abline,
        abline.slope = abline.slope,
        abline.linetype = abline.linetype,
        abline.color = abline.color,
        abline.linewidth = abline.linewidth,
        abline.opacity = abline.opacity,
        legend.show = legend.show,
        legend.color.title = legend.color.title,
        legend.color.breaks = legend.color.breaks,
        legend.color.breaks.labels = legend.color.breaks.labels,
        legend.density.title = legend.density.title,
        legend.density.breaks = legend.density.breaks,
        legend.density.breaks.labels = legend.density.breaks.labels,
        show.grid.lines = show.grid.lines,
        data.out = data.out)
    
    ### RETURN the PLOT ###
    if (data.out) {
        viz_out$to_dittoViz <- pulled_data
        viz_out
    } else {
        viz_out
    }
}
