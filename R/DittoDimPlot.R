
################# dittoDimPlot ####################

#' Shows data overlayed on a tsne, pca, or similar type of plot
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel geom_label_repel
#'
#' @param object A Seurat or SingleCellExperiment object to work with
#' @param var String name of a "gene" or "metadata" (or "ident" for a Seurat \code{object}) to use for coloring the plots.
#' This is the data that will be displayed for each cell/sample.
#' Alternatively, can be a vector of same length as there are cells/samples in the \code{object}.
#' Discrete or continuous data both work. REQUIRED.
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @param reduction.use String, such as "pca", "tsne", "umap", or "PCA", etc, which is the name of a dimensionality reduction slot within the object, and which sets what dimensionality reduction space within the object to use.
#'
#' Default = the first dimensionality reduction slot inside the object named "umap", "tsne", or "pca", or the first dimensionality reduction slot if nonw of those exist.
#' @param size Number which sets the size of data points. Default = 1.
#' @param opacity Number between 0 and 1.
#' Great for when you have MANY overlapping points, this sets how solid the points should be:
#' 1 = not see-through at all. 0 = invisible. Default = 1.
#' (In terms of typical ggplot variables, = alpha)
#' @param dim.1 The component number to use on the x-axis.  Default = 1
#' @param dim.2 The component number to use on the y-axis.  Default = 2
#' @param theme A ggplot theme which will be applied before dittoSeq adjustments. Default = \code{theme_bw()}. See \code{https://ggplot2.tidyverse.org/reference/ggtheme.html} for other options.
#' @param color.panel String vector which sets the colors to draw from. \code{dittoColors()} by default, see \code{\link{dittoColors}} for contents.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param shape.var Variable for setting the shape of cells/samples in the plot.  Note: must be discrete.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#'
#' Note: shapes can be harder to see, and to process mentally, than colors.
#' Even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param shape.panel Vector of integers corresponding to ggplot shapes which sets what shapes to use.
#' When discrete groupings are supplied by \code{shape.var}, this sets the panel of shapes.
#' When nothing is supplied to \code{shape.var}, only the first value is used.
#' Default is a set of 6, \code{c(16,15,17,23,25,8)}, the first being a simple, solid, circle.
#'
#' Note: Unfortunately, shapes can be hard to see when points are on top of each other & they are more slowly processed by the brain.
#' For these reasons, even as a color blind person myself writing this code, I recommend use of colors for variables with many discrete values.
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.size Number representing the size to increase the plotting of color legend shapes to (for discrete variable plotting).
#' Default = 5. *Enlarging the colors legend is incredibly helpful for making colors more distinguishable by color blind individuals.
#' @param legend.title String which sets the title for the main legend which includes the colors. Default = \code{NULL} normally, but \code{var} when a shape legend will also be shown.
#' @param shape.legend.size Number representing the size to increase the plotting of shapes legend shapes to.
#' @param shape.legend.title String which sets the title of the shapes legend.  Default is the \code{shape.var}
#' @param assay,slot single strings or integer that set which data to use when plotting gene expression. See \code{\link{gene}} for more information.
#' @param adjustment When plotting gene expression (or antibody, or other forms of counts data), should that data be used directly (default) or should it be adjusted to be
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
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @param show.others Logical. TRUE by default, whether other cells should be shown in the background in light gray.
#' @param do.ellipse Logical. Whether the groups should be surrounded by an ellipse.
#' @param do.label  Logical. Whether to add text labels at the center (median) of clusters for grouping vars
#' @param labels.size Size of the the labels text
#' @param labels.highlight Logical. Whether the labels should have a box behind them
#' @param labels.repel Logical, that sets whether the labels' placements will be adjusted with \link{ggrepel} to avoid intersections between labels and plot bounds.
#' TRUE by default.
#' @param rename.var.groups String vector which sets new names for the identities of \code{var} groups.
#' @param rename.shape.groups String vector which sets new names for the identities of \code{shape.var} groups.
#' @param min.color color for lowest values of var/min.  Default = yellow
#' @param max.color color for highest values of var/max.  Default = blue
#' @param min Number which sets the value associated with the minimum color.
#' @param max Number which sets the value associated with the maximum color.
#' @param legend.breaks Numeric vector which sets the discrete values to show in the color-scale legend for continuous data.
#' @param legend.breaks.labels String vector, with same length as \code{legend.breaks}, which renames what's displayed next to the tick marks of the color-scale.
#' @param do.letter Logical which sets whether letters should be added on top of the colored dots. For extended colorblindness compatibility.
#' NOTE: \code{do.letter} is ignored if \code{do.hover = TRUE} or \code{shape.var} is provided a metadata because
#' lettering is incompatible with plotly and with changing the dots' to be different shapes.
#' @param do.hover Logical which controls whether the object will be converted to a plotly object so that data about individual points will be displayed when you hover your cursor over them.
#' \code{hover.data} argument is used to determine what data to use.
#' @param hover.data String vector of gene and metadata names, example: \code{c("meta1","gene1","meta2","gene2")} which determines what data to show on hover when \code{do.hover} is set to \code{TRUE}.
#' @param hover.assay,hover.slot,hover.adjustment Similar to the non-hover versions, when showing expression data upon hover, these set what data will be shown.
#' @param add.trajectory.lineages List of vectors representing trajectory paths from start-cluster to end-cluster where vector contents are the names of clusters provided in the \code{trajectory.cluster.meta} input.
#'
#' If the \code{\link[slingshot]{slingshot}} package was used for trajectory analysis, you can use \code{add.trajectory.lineages = SlingshotDataSet(SCE_with_slingshot)$lineages}. In future versions, I might build such retrieval in by default for SCEs.
#' @param add.trajectory.curves List of matrices, each representing coordinates for a trajectory path, from start to end, where matrix columns are x (\code{dim.1}) and y (\code{dim.2}).
#'
#' Alternatively, a list of lists(/princurve objects) can be provided where the target matrices are named as "s" within the list. Thus, if the \code{\link[slingshot]{slingshot}} package was used for trajectory analysis, you can use \code{add.trajectory.curves = SlingshotDataSet(SCE_with_slingshot)$curves}
#' @param trajectory.cluster.meta String name of metadata containing the clusters that were used for generating trajectories.  Required when plotting trajectories using the \code{add.trajectory.lineages} method. Names of clusters inside the metadata should be the same as the contents of \code{add.trajectory.lineages} vectors.
#' @param trajectory.arrow.size Number representing the size of trajectory arrows, in inches.  Default = 0.15.
#' @param data.out Whether just the plot should be output, or a list with the plot and Target_data and Others_data dataframes.  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @return A ggplot or plotly object where colored dots (or other shapes) are overlayed onto a tSNE, PCA, UMAP, ..., plot of choice.
#' Alternatively, a list contatining the ggplot plus the data.frame(s) that went into building it.
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided directly, it just uses that),
#' plus X and Y coordinates data determined by the \code{reduction.use} and \code{dim.1} (x-axis) and \code{dim.2} (y-axis) inputs.
#' The \code{assay}, \code{slot}, and \code{adjustment} inputs can be used to change what expression data is used when displaying gene expression.
#' If a metadata is given to \code{shape.var}, that is retrieved and added to the dataframe as well.
#'
#' Next, if a set of cells or samples to use is indicated with the \code{cells.use} input, then the dataframe is split into \code{Target_data} and \code{Others_data} based on subsetting by the target cells/samples.
#'
#' Finally, a scatter plot is then created using these dataframes where non-target cells will be displayed in gray if \code{show.others=TRUE},
#' and target cell data is displayed on top, colored based on the \code{var}-associated data, and with shapes determined by the \code{shape.var}-associated data.
#'
#' If \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' If \code{do.hover} is set to \code{TRUE}, the plot is coverted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' Many characteristics of the plot can be adjusted using discrete inputs:
#' \itemize{
#' \item \code{size} and \code{opacity} can be used to adjust the size and transparency of the data points.
#' \item Color can be adjusted with \code{color.panel} and/or \code{colors} for discrete data, or \code{min}, \code{max}, \code{min.color}, and \code{max.color} for continuous data.
#' \item Shapes can be adjusted with \code{shape.panel}.
#' \item Color and shape labels can be changed using \code{rename.var.groups} and \code{rename.shape.groups}.
#' \item Titles and axes labels can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item Legends can also be adjusted in other ways, using variables that all start with "\code{legend.}" for easy tab-completion lookup.
#' }
#'
#' Additional Features: Many other tweaks and features can be added as well.
#' Each is accessible through autocompletion starting with "\code{do.}"\code{---} or "\code{add.}"\code{---},
#' and if additional inputs are involved in implementing or tweaking these, the associated inputs will start with the "\code{---.}":
#' \itemize{
#' \item If \code{do.label} is set to \code{TRUE}, labels will be added based on median centers of the discrete \code{var}-data groupings.
#' The size of the text in the labels can be adjusted using the \code{labels.size} input.
#' By default labels will repel eachother and the bounds of the plot, and labels will be highlighted with a white background.
#' Either of these can be turned off by setting \code{labels.repel=FALSE} or \code{labels.highlight=FALSE},
#' \item If \code{do.ellipse} is set to \code{TRUE}, ellipses will be added to highlight distinct \code{var}-data groups' positions based on median positions of their cell/sample components.
#' \item If \code{add.trajectory.lineages} is provided a list of vectors (each vector being cluster names from start-cluster-name to end-cluster-name), and a metadata name pointing to the relevant clustering information is provided to \code{trajectory.cluster.meta},
#' then median centers of the clusters will be calculated and arrows will be overlayed to show trajectory inference paths in the current dimmenionality reduction space.
#' \item If \code{add.trajectory.curves} is provided a list of matrices (each matrix containing x, y coordinates from start to end), paths and arrows will be overlayed to show trajectory inference curves in the current dimmenionality reduction space.
#' Arrow size is controlled with the \code{trajectory.arrow.size} input.
#' }
#'
#' @seealso
#' \code{\link{getGenes}} and \code{\link{getMetas}} to see what the \code{var}, \code{shape.var}, and \code{hover.data} options are.
#'
#' \code{\link{importDittoBulk}} for how to create a \code{\link{SingleCellExperiment}} object from bulk seq data that dittoSeq functions can use &
#' \code{\link{addDimReduction}} for how to add calculated dimensionality reductions that \code{dittoDimPlot} can utilize.
#'
#' \code{\link{dittoScatterPlot}} for showing very similar data representations, but where genes or metadata are wanted as the axes.
#'
#' \code{\link{dittoPlot}} for an alternative continuous data display method where data is shown on a y- (or x-) axis.
#'
#' \code{\link{dittoBarPlot}} for an alternative discrete data display and quantification method.
#'
#' @author Daniel Bunis
#' @importFrom stats median
#' @export
#' @examples
#' pbmc <- Seurat::pbmc_small
#'
#' # Display discrete data:
#' dittoDimPlot(pbmc, "RNA_snn_res.1")
#' # Display continuous data:
#' dittoDimPlot(pbmc, "CD14")
#'
#' # To show currently set clustering for seurat objects, you can use "ident".
#' # To change the dimensional reduction type, use reduction.use.
#'
#' # MANY addtional tweaks are possible.
#' # Also, many extra features are easy to add as well:
#' dittoDimPlot(pbmc, "ident", do.label = TRUE)
#' dittoDimPlot(pbmc, "ident", do.label = TRUE, do.ellipse = TRUE)
#' dittoDimPlot(pbmc, "CD3E", do.hover = TRUE,
#'     hover.data = c("CD14", "RNA_snn_res.0.8", "groups"))
#' dittoDimPlot(pbmc, "CD3E", add.trajectory.lineages = list(c(0:2), c(0,2)),
#'     trajectory.cluster.meta = "ident")

dittoDimPlot <- function(
    object, var="ident", reduction.use = NA, size=1, opacity = 1,
    dim.1 = 1, dim.2 = 2, cells.use = NULL, show.others=TRUE,
    show.axes.numbers = TRUE,
    color.panel = dittoColors(), colors = seq_along(color.panel),
    shape.var = NULL, shape.panel=c(16,15,17,23,25,8),
    assay = .default_assay(object), slot = .default_slot(object),
    adjustment = NULL,
    main = "make", sub = NULL, xlab = "make", ylab = "make",
    theme = NA, legend.show = TRUE, legend.size = 5,
    legend.title = "make",
    shape.legend.size = 5, shape.legend.title = shape.var,
    do.ellipse = FALSE, do.label = FALSE,
    labels.size = 5, labels.highlight = TRUE, labels.repel = TRUE,
    rename.var.groups = NULL, rename.shape.groups = NULL,
    min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
    legend.breaks = waiver(), legend.breaks.labels = waiver(),
    do.letter = FALSE, do.hover = FALSE, hover.data = var,
    hover.assay = .default_assay(object), hover.slot = .default_slot(object),
    hover.adjustment = NULL,
    add.trajectory.lineages = NULL, add.trajectory.curves = NULL,
    trajectory.cluster.meta, trajectory.arrow.size = 0.15, data.out = FALSE) {

    #Standardize cells.use to a list of names.
    cells.use <- .which_cells(cells.use, object)
    all.cells <- .all_cells(object)

    if (do.hover || !is.null(shape.var)) {
        do.letter <- FALSE
    }

    # Generate the x/y dimensional reduction data and plot titles.
    if (is.na(reduction.use)) {
        reduction.use <- .default_reduction(object)
    }
    xdat <- .extract_Reduced_Dim(reduction.use, dim.1, object)
    ydat <- .extract_Reduced_Dim(reduction.use, dim.2, object)
    xlab <- .leave_default_or_null(xlab, xdat$name)
    ylab <- .leave_default_or_null(ylab, ydat$name)
    main <- .leave_default_or_null(main, var, length(var)!=1)
    legend.title <- .leave_default_or_null(legend.title, var, is.null(shape.var))

    # Edit theme.
    if (is.na(theme[1])){
        theme <- theme_bw()
        if (grepl("tsne|umap", tolower(reduction.use))) {
            # Remove grid lines
            theme <- theme +
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())
        }
    }
    if (!show.axes.numbers) {
        theme <- theme +
            theme(axis.text.x=element_blank(),axis.text.y=element_blank())
    }

    # Make dataframes and plot
    p.df <- dittoScatterPlot(
        xdat$embeddings, ydat$embeddings, var, shape.var, object, cells.use,
        show.others, size, opacity, color.panel, colors,
        NULL, NULL, NULL, NULL, NULL, NULL, assay, slot, adjustment,
        do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
        shape.panel, rename.var.groups, rename.shape.groups,
        min.color, max.color, min, max,
        xlab, ylab, main, sub, theme, legend.show, legend.title, legend.size,
        legend.breaks, legend.breaks.labels, shape.legend.title,
        shape.legend.size, data.out = TRUE)
    p <- p.df$plot
    Target_data <- p.df$Target_data
    Others_data <- p.df$Others_data

    # Add extra features
    if (do.letter) {
        p <- .add_letters(
            p, Target_data, "color", size, opacity, legend.title, legend.size)
    }
    if (do.ellipse) {
        p <- p + stat_ellipse(
            data=Target_data,
            aes_string(x = "X", y = "Y", colour = "color"),
            type = "t", linetype = 2, size = 0.5, show.legend = FALSE)
    }
    if (do.label) {
        p <- .add_labels(
            p, Target_data, "color", labels.highlight, labels.size,
            labels.repel)
    }
    if (is.list(add.trajectory.lineages)) {
        p <- .add_trajectory_lineages(
            p, add.trajectory.lineages, trajectory.cluster.meta,
            trajectory.arrow.size, object, reduction.use, dim.1, dim.2)
    }

    if (is.list(add.trajectory.curves)) {
        p <- .add_trajectory_curves(
            p, add.trajectory.curves, trajectory.arrow.size, dim.1, dim.2)
    }

    if (!legend.show) {
        p <- .remove_legend(p)
    }
    ### RETURN the PLOT ###
    if (data.out) {
        return(list(
            plot = p,
            Target_data = Target_data,
            Others_data = Others_data))
    } else {
        if (do.hover) {
            .error_if_no_plotly()
            return(plotly::ggplotly(p, tooltip = "text"))
        } else {
            return(p)
        }
    }
}

.default_reduction <- function(object) {
    # Use umap > tsne > pca, or whatever the first reduction slot is.
    opts <- getReductions(object)
    if (is.null(opts)) {
        stop("No dimensionality reductions available.")
    }
    use <- .preferred_or_first(opts, c("umap","tsne","pca"))
    use
}

.add_labels <- function(p, Target_data, col.use = "color", labels.highlight, labels.size, labels.repel) {
    #Make a text plot at the median x and y values for each cluster

    #Determine medians
    cent.x = vapply(
        levels(as.factor(Target_data[,col.use])),
        function(level) {
            median(Target_data$X[Target_data[,col.use]==level])
        }, FUN.VALUE = numeric(1))
    cent.y = vapply(
        levels(as.factor(Target_data[,col.use])),
        function(level) {
            median(Target_data$Y[Target_data[,col.use]==level])
        }, FUN.VALUE = numeric(1))

    #Add labels
    args <- list(
        data = data.frame(cent.x=cent.x, cent.y=cent.y),
        mapping = aes(x = cent.x, y = cent.y),
        size = labels.size,
        label = levels(as.factor(Target_data[,col.use])))
    geom.use <-
        if (labels.highlight) {
            if (labels.repel) {
                ggrepel::geom_label_repel
            } else {
                geom_label
            }
        } else {
            if (labels.repel) {
                ggrepel::geom_text_repel
            } else {
                geom_text
            }
        }
    p + do.call(geom.use, args)
}

.add_trajectory_lineages <- function(
    p, trajectories, clusters, arrow.size = 0.15, object, reduction.use,
    dim.1, dim.2) {
    # p = the $p output of a dittoDimPlot(any.var,..., data.out = TRUE)
    # clusters = the name of the metadata metadata slot that holds the clusters used for cluster-based trajectory analysis
    # trajectories = List of lists of cluster paths. Also, the output of SlingshotDataSet(SCE_with_slingshot)$lineages
    # arrow.size = numeric scalar that sets the arrow length (in inches) at the endpoints of trajectory lines.

    #Determine medians
    cluster.dat <- meta(clusters, object)
    cluster.levels <- meta.levels(clusters, object)
    data <- data.frame(
        cent.x = vapply(
            cluster.levels,
            function(level) {
                median(
                    .extract_Reduced_Dim(reduction.use, dim.1, object)$embedding[
                        cluster.dat==level])
            }, FUN.VALUE = numeric(1)),
        cent.y = vapply(
            cluster.levels,
            function(level) {
                median(
                    .extract_Reduced_Dim(reduction.use, dim.2, object)$embedding[
                        cluster.dat==level])
            }, FUN.VALUE = numeric(1)))

    #Add trajectories
    for (i in seq_along(trajectories)){
        p <- p + geom_path(
            data = data[as.character(trajectories[[i]]),],
            aes_string(x = "cent.x", y = "cent.y"),
            arrow = arrow(
                angle = 20, type = "closed", length = unit(arrow.size, "inches")))
    }
    p
}

.add_trajectory_curves <- function(
    p, trajectories, arrow.size = 0.15, dim.1, dim.2) {
    # p = the $p output of a dittoDimPlot(any.var,..., data.out = TRUE)
    # trajectories = List of matrices containing trajectory curves. The output of SlingshotDataSet(SCE_with_slingshot)$curves can be used if the coordinate matrix (`$s`) for each list is extracted and they are all stored in a list.
    # arrow.size = numeric scalar that sets the arrow length (in inches) at the endpoints of trajectory lines.

    if ("s" %in% names(trajectories[[1]])) {
    #Add trajectories for princurves/slingshot list of lists provision method
        for (i in seq_along(trajectories)){
            #extract fit coords per cell
            data <- as.data.frame(trajectories[[i]]$s)
            #order cells' fit coords by pseudotime order
            data <- data[trajectories[[i]]$ord,]
            #name the dimensions used
            names(data)[c(dim.1,dim.2)] <- c("x", "y")
            p <- p + geom_path(
                data = data,
                aes_string(x = "x", y = "y"),
                arrow = arrow(
                    angle = 20, type = "closed", length = unit(arrow.size, "inches")))
        }
    } else {
    #Add trajectories for general list of matrices provision method.
    #  Note: Accepts dataframes too.
        for (i in seq_along(trajectories)){
            data <- as.data.frame(trajectories[[i]])
            names(data) <- c("x", "y")
            p <- p + geom_path(
                data = data,
                aes_string(x = "x", y = "y"),
                arrow = arrow(
                    angle = 20, type = "closed", length = unit(arrow.size, "inches")))
        }
    }
    p
}

.add_letters <- function(
    p, Target_data, col.use = "color", size, opacity, legend.title,
    legend.size) {

    letters.needed <- length(levels(as.factor(Target_data[,col.use])))
    letter.labels <- c(
        LETTERS, letters, 0:9, "!", "@", "#", "$", "%", "^", "&", "*", "(",
        ")", "-", "+", "_", "=", ";", "/", "|", "{", "}", "~"
        )[seq_len(letters.needed)]
    names(letter.labels) <- levels(as.factor(Target_data[,col.use]))
    p <- p +
        geom_point(
            data=Target_data,
            aes_string(x = "X", y = "Y", shape = col.use),
            color = "black", size=size*3/4, alpha = opacity) +
        scale_shape_manual(
            name = legend.title,
            values = letter.labels)
    p
}

#### multi_dittoDimPlot : a function for quickly making multiple DBDimPlots arranged in a grid.
#' Generates multiple dittoDimPlots arranged in a grid.
#'
#' @param object A Seurat or SingleCellExperiment object to work with
#' @param vars c("var1","var2","var3",...). A list of vars from which to generate the separate plots
#' @param ncol,nrow Integer/NULL. How many columns or rows the plots should be arranged into
#' @param axes.labels.show Logical. Whether a axis labels should be shown. Ignored if xlab or ylab are set manually.
#' @param OUT.List Logical. (Default = FALSE) When set to \code{TRUE}, a list of the individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' @param legend.show,xlab,ylab,... other paramters passed to dittoDimPlot.
#' @return Given multiple 'var' parameters to \code{vars}, this function will output a dittoDimPlot for each one, arranged into a grid, with some slight tweaks to the defaults.
#' If \code{OUT.list} was set to TRUE, the list of individual plots, named by the \code{vars} being shown in each, is output instead of the combined multi-plot.
#' All parameters that can be adjusted in dittoDimPlot can be adjusted here, but the only parameter that can be adjusted between each is the \code{var}.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#'
#' genes <- c("CD8A","CD3E","FCER1A","CD14","MS4A1")
#' multi_dittoDimPlot(c(genes, "ident"), object = "pbmc")
#'
#' @author Daniel Bunis
#' @export

multi_dittoDimPlot <- function(
    object, vars, legend.show = FALSE, ncol = NULL, nrow = NULL,
    axes.labels.show = FALSE, xlab = NA, ylab = NA, OUT.List = FALSE, ...) {

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
