
################# dittoDimPlot ####################

#' Shows data overlayed on a tsne, pca, or similar type of plot
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel geom_label_repel
#'
#' @param object A Seurat, SingleCellExperiment, or SummarizedExperiment object.
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample. Discrete or continuous data both work.
#'
#' Alternatively, a string vector naming multiple genes or metadata, OR a vector of the same length as there are cells/samples in the \code{object} which provides per-cell data directly.
#' @param reduction.use String, such as "pca", "tsne", "umap", or "PCA", etc, which is the name of a dimensionality reduction slot within the object, and which sets what dimensionality reduction space within the object to use.
#'
#' Default = the first dimensionality reduction slot inside the object with "umap", "tsne", or "pca" within its name, (priority: UMAP > t-SNE > PCA) or the first dimensionality reduction slot if none of those exist.
#' 
#' Alternatively, a matrix (or data.frame) containing the dimensionality reduction embeddings themselves.
#' The matrix should have as many rows as there are cells/samples in the \code{object}.
#' Note that \code{dim.1} and \code{dim.2} will still be used to select which columns to pull from, and column names will serve as the default \code{xlab} & \code{ylab}.
#' @param size Number which sets the size of data points. Default = 1.
#' @param opacity Number between 0 and 1.
#' Great for when you have MANY overlapping points, this sets how solid the points should be:
#' 1 = not see-through at all. 0 = invisible. Default = 1.
#' (In terms of typical ggplot variables, = alpha)
#' @param dim.1 The component number to use on the x-axis. Default = 1
#' @param dim.2 The component number to use on the y-axis. Default = 2
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments.
#' Default = \code{theme_bw()}.
#' See \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options and ideas.
#' @param show.grid.lines Logical which sets whether gridlines of the plot should be shown.
#' They are removed when set to FALSE.
#' Default = FALSE for umap and tsne \code{reduction.use}, TRUE otherwise.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default, see \code{\link{dittoColors}} for contents.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use.
#'
#' Useful for quickly swapping the colors of nearby clusters.
#' @param shape.by Variable for setting the shape of cells/samples in the plot.  Note: must be discrete.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#'
#' Note: shapes can be harder to see, and to process mentally, than colors.
#' Even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param multivar.split.dir "row" or "col", sets the direction of faceting used for 'var' values when \code{var} is given multiple genes or metadata, and when \code{split.by} is used to provide additional data to facet by.
#' @param split.by 1 or 2 strings naming discrete metadata to use for splitting the cells/samples into multiple plots with ggplot faceting.
#'
#' When 2 metadatas are named, c(row,col), the first is used as rows and the second is used for columns of the resulting grid.
#'
#' When 1 metadata is named, shape control can be achieved with \code{split.nrow} and \code{split.ncol}
#'
#' @param split.nrow,split.ncol Integers which set the dimensions of faceting/splitting when a single metadata is given to \code{split.by}.
#' @param split.show.all.others Logical which sets whether gray "others" cells of facets should include all cells of other facets (\code{TRUE}) versus just cells left out by \code{cell.use} (\code{FALSE}).
#' @param split.adjust A named list which allows extra parameters to be pushed through to the faceting function call.
#' List elements should be valid inputs to the faceting functions, e.g. `list(scales = "free")`.
#' 
#' For options, when giving 1 metadata to \code{split.by}, see \code{\link[ggplot2]{facet_wrap}},
#' OR when giving 2 metadatas to \code{split.by}, see \code{\link[ggplot2]{facet_grid}}.
#' @param extra.vars String vector providing names of any extra metadata to be stashed in the dataframe supplied to \code{ggplot(data)}.
#'
#' Useful for making custom splitting/faceting or other additional alterations \emph{after} dittoSeq plot generation.
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.by}, this sets the panel of shapes.
#' When nothing is supplied to \code{shape.by}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#'
#' Note: Unfortunately, shapes can be hard to see when points are on top of each other & they are more slowly processed by the brain.
#' For these reasons, even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.size Number representing the size at which color legend shapes should be plotted (for discrete variable plotting) in the color legend.
#' Default = 5. *Enlarging the colors legend is incredibly helpful for making colors more distinguishable by color blind individuals.
#' @param legend.title String which sets the title for the color legend. Default = \code{NULL} normally, but \code{var} when a shape legend will also be shown.
#' @param shape.legend.size Number representing the size at which shapes should be plotted in the shape legend.
#' @param shape.legend.title String which sets the title of the shapes legend.  Default is \code{shape.by}
#' @param assay,slot single strings or integers (SCEs and SEs) or an optionally named vector of such values that set which expression data to use.
#' See \code{\link{GeneTargeting}} for specifics and examples -- Seurat and SingleCellExperiment objects deal with these differently, and functionality additions in dittoSeq have led to some minimal divergence from the native methodologies.
#' @param adjustment When plotting gene / feature expression, should that data be used directly (default) or should it be adjusted to be
#' \itemize{
#' \item{"z-score": scaled with the scale() function to produce a relative-to-mean z-score representation}
#' \item{"relative.to.max": divided by the maximum expression value to give percent of max values between [0,1]}
#' }
#' @param main String, sets the plot title.
#' Default title is automatically generated if not given a specific value.  To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle
#' @param xlab,ylab Strings which set the labels for the axes.
#' Default labels are generated if you do not give this a specific value.
#' To remove, set to \code{NULL}.
#' @param show.axes.numbers Logical which controls whether the axes values should be displayed.
#' @param cells.use String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' @param show.others Logical. Whether other cells should be shown in the background in light gray. Default = TRUE.
#' @param do.ellipse Logical. Whether the groups should be surrounded by median-centered ellipses.
#' @param do.label  Logical. Whether to add text labels near the center (median) of clusters for grouping vars.
#' @param labels.size Size of the the labels text
#' @param labels.highlight Logical. Whether the labels should have a box behind them
#' @param labels.repel Logical, that sets whether the labels' placements will be adjusted with \link{ggrepel} to avoid intersections between labels and plot bounds.
#' TRUE by default.
#' @param labels.split.by String of one or two metadata names which controls the facet-split calculations for label placements.
#' Defaults to \code{split.by}, so generally there is no need to adjust this except when you are utilizing the \code{extra.vars} input to achieve manual faceting control.
#' @param labels.repel.adjust A named list which allows extra parameters to be pushed through to ggrepel function calls.
#' List elements should be valid inputs to the \code{\link[ggrepel]{geom_label_repel}} by default, or \code{\link[ggrepel]{geom_text_repel}} when \code{labels.highlight = FALSE}.
#' @param rename.var.groups String vector which sets new names for the identities of \code{var} groups.
#' @param rename.shape.groups String vector which sets new names for the identities of \code{shape.by} groups.
#' @param min.color color for lowest values of \code{var}/\code{min}.  Default = yellow
#' @param max.color color for highest values of \code{var}/\code{max}.  Default = blue
#' @param min Number which sets the value associated with the minimum color.
#' @param max Number which sets the value associated with the maximum color.
#' @param order String. If the data should be plotted based on the order of the color data, sets whether to plot (from back to front) in "increasing", "decreasing", "randomize" order.
#' If left as "unordered", plot order is simply based on the order of cells within the \code{object}.
#' @param legend.breaks Numeric vector which sets the discrete values to show in the color-scale legend for continuous data.
#' @param legend.breaks.labels String vector, with same length as \code{legend.breaks}, which renames what's displayed next to the tick marks of the color-scale.
#' @param do.letter Logical which sets whether letters should be added on top of the colored dots. For extended colorblindness compatibility.
#' NOTE: \code{do.letter} is ignored if \code{do.hover = TRUE} or \code{shape.by} is provided a metadata because
#' lettering is incompatible with plotly and with changing the dots' to be different shapes.
#' @param do.contour Logical. Whether density-based contours should be displayed.
#' @param contour.color String that sets the color(s) of the \code{do.contour} contours.
#' @param contour.linetype String or numeric which sets the type of line used for \code{do.contour} contours.
#' Defaults to "solid", but see \code{\link[ggplot2]{linetype}} for other options.
#' @param do.hover Logical which controls whether the output will be converted to a plotly object so that data about individual points will be displayed when you hover your cursor over them.
#' \code{hover.data} argument is used to determine what data to use.
#' @param hover.data String vector of gene and metadata names, example: \code{c("meta1","gene1","meta2")} which determines what data to show on hover when \code{do.hover} is set to \code{TRUE}.
#' @param hover.assay,hover.slot,hover.adjustment Similar to the non-hover versions of these inputs, when showing expression data upon hover, these set what data will be shown.
#' @param add.trajectory.lineages List of vectors representing trajectory paths, each from start-cluster to end-cluster, where vector contents are the names of clusters provided in the \code{trajectory.cluster.meta} input.
#'
#' If the \code{\link[slingshot]{slingshot}} package was used for trajectory analysis,
#' you can provide \code{add.trajectory.lineages = slingLineages('object')}.
#' @param add.trajectory.curves List of matrices, each representing coordinates for a trajectory path, from start to end, where matrix columns represent x (\code{dim.1}) and y (\code{dim.2}) coordinates of the paths.
#'
#' Alternatively, a list of lists(/princurve objects) can be provided.
#' Thus, if the \code{\link[slingshot]{slingshot}} package was used for trajectory analysis,
#' you can provide \code{add.trajectory.curves = slingCurves('object')}
#' @param trajectory.cluster.meta String name of metadata containing the clusters that were used for generating trajectories.  Required when plotting trajectories using the \code{add.trajectory.lineages} method. Names of clusters inside the metadata should be the same as the contents of \code{add.trajectory.lineages} vectors.
#' @param trajectory.arrow.size Number representing the size of trajectory arrows, in inches.  Default = 0.15.
#' @param do.raster Logical. When set to \code{TRUE}, rasterizes the internal plot area. Useful for editing in external programs (e.g. Illustrator).
#' @param raster.dpi Number indicating dpi to use for rasterization. Default = 300.
#' @param data.out Logical. When set to \code{TRUE}, changes the output, from the plot alone, to a list containing the plot ("p"),
#' a data.frame containing the underlying data for target cells ("Target_data"),
#' and a data.frame containing the underlying data for non-target cells ("Others_data").
#'
#' @inheritParams gene
#'
#' @return A ggplot or plotly object where colored dots (or other shapes) are overlayed onto a tSNE, PCA, UMAP, ..., plot of choice.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' Alternatively, if \code{do.hover} is set to \code{TRUE}, the plot is coverted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided directly, it just uses that),
#' plus X and Y coordinates data determined by the \code{reduction.use} and \code{dim.1} (x-axis) and \code{dim.2} (y-axis) inputs.
#' Any extra data requested with \code{shape.by}, \code{split.by} or \code{extra.var} is added as well.
#' For expression/counts data, \code{assay}, \code{slot}, and \code{adjustment} inputs can be used to change which data is used, and if it should be adjusted in some way.
#'
#' Next, if a set of cells or samples to use is indicated with the \code{cells.use} input, then the dataframe is split into \code{Target_data} and \code{Others_data} based on subsetting by the target cells/samples.
#'
#' Finally, a scatter plot is then created using these dataframes where non-target cells will be displayed in gray if \code{show.others=TRUE},
#' and target cell data is displayed on top, colored based on the \code{var}-associated data, and with shapes determined by the \code{shape.by}-associated data.
#' If \code{split.by} was used, the plot will be split into a matrix of panels based on the associated groupings.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item \code{size} and \code{opacity} can be used to adjust the size and transparency of the data points.
#' \item Color can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes can be adjusted with \code{shape.panel}.
#' \item Color and shape labels can be changed using \code{rename.var.groups} and \code{rename.shape.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab-completion lookup.
#' }
#'
#' @section Additional Features:
#' Many other tweaks and features can be added as well.
#' Each is accessible through 'tab' autocompletion starting with "\code{do.}"\code{---} or "\code{add.}"\code{---},
#' and if additional inputs are involved in implementing or tweaking these, the associated inputs will start with the "\code{---.}":
#' \itemize{
#' \item If \code{do.label} is set to \code{TRUE}, labels will be added based on median centers of the discrete \code{var}-data groupings.
#' The size of the text in the labels can be adjusted using the \code{labels.size} input.
#' By default labels will repel eachother and the bounds of the plot, and labels will be highlighted with a white background.
#' Either of these can be turned off by setting \code{labels.repel = FALSE} or \code{labels.highlight = FALSE},
#' \item If \code{do.ellipse} is set to \code{TRUE}, ellipses will be added to highlight distinct \code{var}-data groups' positions based on median positions of their cell/sample components.
#' \item If \code{do.contour} is provided, density gradiant contour lines will be overlaid with color and linetype adjustable via \code{contour.color} and \code{contour.linetype}.
#' \item If \code{add.trajectory.lineages} is provided a list of vectors (each vector being cluster names from start-cluster-name to end-cluster-name), and a metadata name pointing to the relevant clustering information is provided to \code{trajectory.cluster.meta},
#' then median centers of the clusters will be calculated and arrows will be overlayed to show trajectory inference paths in the current dimmenionality reduction space.
#' \item If \code{add.trajectory.curves} is provided a list of matrices (each matrix containing x, y coordinates from start to end), paths and arrows will be overlayed to show trajectory inference curves in the current dimmenionality reduction space.
#' Arrow size is controlled with the \code{trajectory.arrow.size} input.
#' }
#'
#' @seealso
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{var}, \code{split.by}, etc. options are of an \code{object}.
#'
#' \code{\link{getReductions}} to see what the \code{reduction.use} options are of an \code{object}.
#' 
#' \code{\link{importDittoBulk}} for how to create a \code{\link{SingleCellExperiment}} object from bulk seq data that dittoSeq functions can use &
#' \code{\link{addDimReduction}} for how to specifically add calculated dimensionality reductions that \code{dittoDimPlot} can utilize.
#'
#' \code{\link{dittoScatterPlot}} for showing very similar data representations, but where genes or metadata are wanted as the axes.
#' 
#' \code{\link{dittoDimHex}} and \code{\link{dittoScatterHex}} for showing very similar data representations, but where nearby cells are summarized together in hexagonal bins.
#'
#' \code{\link{dittoPlot}} for an alternative continuous data display method where data broken into discrete groupings is shown on a y- (or x-) axis.
#'
#' \code{\link{dittoBarPlot}} for an alternative discrete data display and quantification method.
#'
#' @author Daniel Bunis and Jared Andrews
#' @importFrom stats median
#' @export
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' # Display discrete data:
#' dittoDimPlot(myRNA, "clustering")
#' # Display continuous data:
#' dittoDimPlot(myRNA, "gene1")
#' 
#' # You can also plot multiple sets of continuous data:
#' dittoDimPlot(myRNA, c("gene1", "gene2"))
#' # (See ?multi_dittoDimPlot if you would like to have wholy separate
#' # plots/scales/legends for each set.)
#'
#' # To show currently set clustering for seurat objects, you can use "ident".
#' # To change the dimensional reduction type, use 'reduction.use'.
#' dittoDimPlot(myRNA, "clustering",
#'     reduction.use = "pca",
#'     dim.1 = 3,
#'     dim.2 = 4)
#'
#' # Subset to certain cells with cells.use
#' dittoDimPlot(myRNA, "clustering",
#'     cells.us = !myRNA$SNP)
#'
#' # Data can also be split in other ways with 'shape.by' or 'split.by'
#' dittoDimPlot(myRNA, "gene1",
#'     shape.by = "clustering",
#'     split.by = "SNP") # single split.by element
#' dittoDimPlot(myRNA, "gene1",
#'     split.by = c("groups","SNP")) # row and col split.by elements
#'
#' # Modify the look with intuitive inputs
#' dittoDimPlot(myRNA, "clustering",
#'     size = 2, opacity = 0.7, show.axes.numbers = FALSE,
#'     ylab = NULL, xlab = "tSNE",
#'     main = "Plot Title",
#'     sub = "subtitle",
#'     legend.title = "clustering")
#'
#' # MANY addtional tweaks are possible.
#' # Also, many extra features are easy to add as well:
#' dittoDimPlot(myRNA, "clustering",
#'     do.label = TRUE, do.ellipse = TRUE)
#' dittoDimPlot(myRNA, "clustering",
#'     do.label = TRUE, labels.highlight = FALSE, labels.size = 8)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'     dittoDimPlot(myRNA, "gene1", do.hover = TRUE,
#'         hover.data = c("gene2", "clustering", "timepoint"))
#' }
#' dittoDimPlot(myRNA, "gene1", add.trajectory.lineages = list(c(1,2,4), c(1,3)),
#'     trajectory.cluster.meta = "clustering",
#'     sub = "Pseudotime Trajectories")
#' 
#' dittoDimPlot(myRNA, "gene1",
#'     do.contour = TRUE,
#'     contour.color = "lightblue", # Optional, black by default
#'     contour.linetype = "dashed") # Optional, solid by default
#' 
#' # Plotting ordering can also be adjusted with 'order':
#' dittoDimPlot(myRNA, "timepoint", size = 20,
#'     order = "increasing")
#' dittoDimPlot(myRNA, "timepoint", size = 20,
#'     order = "decreasing")
#' dittoDimPlot(myRNA, "timepoint", size = 20,
#'     order = "randomize")

dittoDimPlot <- function(
    object,
    var,
    reduction.use = .default_reduction(object),
    size = 1,
    opacity = 1,
    dim.1 = 1,
    dim.2 = 2,
    cells.use = NULL,
    shape.by = NULL,
    split.by = NULL,
    split.adjust = list(),
    extra.vars = NULL,
    multivar.split.dir = c("col", "row"),
    show.others = TRUE,
    split.show.all.others = TRUE,
    split.nrow = NULL,
    split.ncol = NULL,
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    shape.panel = c(16,15,17,23,25,8),
    min.color = "#F0E442",
    max.color = "#0072B2",
    min = NA,
    max = NA,
    order = c("unordered", "increasing", "decreasing", "randomize"),
    main = "make",
    sub = NULL,
    xlab = "make",
    ylab = "make",
    rename.var.groups = NULL,
    rename.shape.groups = NULL,
    theme = theme_bw(),
    show.axes.numbers = TRUE,
    show.grid.lines = if (is.character(reduction.use)) { !grepl("umap|tsne", tolower(reduction.use)) } else {TRUE},
    do.letter = FALSE,
    do.ellipse = FALSE,
    do.label = FALSE,
    labels.size = 5,
    labels.highlight = TRUE,
    labels.repel = TRUE,
    labels.split.by = split.by,
    labels.repel.adjust = list(),
    do.hover = FALSE,
    hover.data = var,
    hover.assay = .default_assay(object),
    hover.slot = .default_slot(object),
    hover.adjustment = NULL,
    add.trajectory.lineages = NULL,
    add.trajectory.curves = NULL,
    trajectory.cluster.meta,
    trajectory.arrow.size = 0.15,
    do.contour = FALSE,
    contour.color = "black",
    contour.linetype = 1,
    legend.show = TRUE,
    legend.size = 5,
    legend.title = "make",
    legend.breaks = waiver(),
    legend.breaks.labels = waiver(),
    shape.legend.size = 5,
    shape.legend.title = shape.by,
    do.raster = FALSE,
    raster.dpi = 300,
    data.out = FALSE) {

    order <- match.arg(order)
    multivar.split.dir <- match.arg(multivar.split.dir)
    
    if (do.hover || !is.null(shape.by)) {
        do.letter <- FALSE
    }

    # Generate the x/y dimensional reduction data and plot titles.
    xdat <- .extract_Reduced_Dim(reduction.use, dim.1, object)
    ydat <- .extract_Reduced_Dim(reduction.use, dim.2, object)
    xlab <- .leave_default_or_null(xlab, xdat$name)
    ylab <- .leave_default_or_null(ylab, ydat$name)
    main <- .leave_default_or_null(main, var, length(var)>1)
    legend.title <- .leave_default_or_null(
        legend.title, var, is.null(shape.by) || length(var)>1)

    # Edit theme.
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
    p.df <- dittoScatterPlot(
        object, xdat$embeddings, ydat$embeddings, var, shape.by, split.by,
        extra.vars, cells.use, multivar.split.dir, show.others, split.show.all.others,
        size, opacity, color.panel, colors,
        split.nrow, split.ncol, split.adjust, NA, NA, NA, NA, NA, NA,
        assay, slot, adjustment, assay, slot, adjustment, swap.rownames,
        shape.panel, rename.var.groups, rename.shape.groups,
        min.color, max.color, min, max, order,
        xlab, ylab, main, sub, theme,
        do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
        do.contour, contour.color, contour.linetype,
        add.trajectory.lineages, add.trajectory.curves = NULL,
        trajectory.cluster.meta, trajectory.arrow.size,
        do.letter, do.ellipse, do.label, labels.size, labels.highlight,
        labels.repel, labels.split.by, labels.repel.adjust,
        legend.show, legend.title, legend.size,
        legend.breaks, legend.breaks.labels, shape.legend.title,
        shape.legend.size, do.raster, raster.dpi, data.out = TRUE)
    p <- p.df$plot
    Target_data <- p.df$Target_data
    Others_data <- p.df$Others_data

    # Add extra features
    if (is.list(add.trajectory.curves)) {
        p <- .add_trajectory_curves(
            p, add.trajectory.curves, trajectory.arrow.size, dim.1, dim.2)
    }

    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    if (do.hover) {
        .error_if_no_plotly()
        p <- plotly::ggplotly(p, tooltip = "text")
    }
    
    ### RETURN the PLOT ###
    if (data.out) {
        list(
            plot = p,
            Target_data = Target_data,
            Others_data = Others_data)
    } else {
        p
    }
}
