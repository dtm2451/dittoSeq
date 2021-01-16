#' Plot spatial transcriptomics data by spatial coordinates
#' @inheritParams dittoDimPlot
#' @param split.by ...
#' @param samples.use ...
#' @param images.use ...
#' @return a ggplot scatterplot where colored dots and/or shapes represent data from the given spatial coordinates.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' Alternatively, if \code{do.hover} is set to \code{TRUE} and no images are to be shown (incompatibility), the plot is converted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' @details X/Y are extracted...
#' @export
#' @author Daniel Bunis


dittoSpatial <- function(
    object,
    var = NULL,
    split.by = "sample_id",
    samples.use = NULL,
    images.use = if (is(object, "VisiumExperiment")) {
        imagePaths(object)[1] } else { NULL },
    cells.use = NULL,
    size = 1,
    opacity = 1,
    shape.by = NULL,
    extra.vars = NULL,
    split.nrow = NULL,
    split.ncol = NULL,
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    shape.panel = c(16,15,17,23,25,8),
    show.others = TRUE,
    show.axes.numbers = TRUE,
    show.grid.lines = if (is.character(reduction.use)) { !grepl("umap|tsne", tolower(reduction.use)) } else {TRUE},
    min.color = "#F0E442",
    max.color = "#0072B2",
    min = NULL,
    max = NULL,
    order = c("unordered", "increasing", "decreasing"),
    main = "make",
    sub = NULL,
    xlab = NULL,
    ylab = NULL,
    rename.var.groups = NULL,
    rename.shape.groups = NULL,
    theme = theme_bw(),
    do.letter = FALSE,
    do.ellipse = FALSE,
    do.label = FALSE,
    labels.size = 5,
    labels.highlight = TRUE,
    labels.repel = TRUE,
    labels.split.by = split.by,
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
    data.out = FALSE
    ) {

    is_visium <- is(object,"VisiumExperiment")
    
    order <- match.arg(order)
    
    if (do.hover || !is.null(shape.by)) {
        do.letter <- FALSE
    }

    # Generate the x/y dimensional reduction data and plot titles.
    main <- .leave_default_or_null(main, var, length(var)!=1)
    legend.title <- .leave_default_or_null(legend.title, var, is.null(shape.by))

    # Retrieve coordinates data
    coords <- SpatialExperiment::spatialCoords(object)
    
    if (is_visium) {
        xdat <- coords$array_row*scaleFactors(object)
        ydat <- coords$array_col*scaleFactors(object)
    } else {
        xdat <- coords$x
        ydat <- coords$y
    }

    # Retrieve image data
    image_paths <- images.use
    
    # Determine relevant images per sample if given as a type.
    # if images.use %in% c()
    
    images <- lapply(image_paths, function(x) {
        readbitmap::read.bitmap(x)
        })
    
    images_tibble <- dplyr::tibble(
        height = lapply(images, nrow),
        width = lapply(images, ncol),
        grob = lapply(images, function(x) {
            grid::rasterGrob(x, width=unit(1,"npc"), height=unit(1,"npc"))
        })
    )
    
    ### Subset by sample_id
    # cells.use
    # images_tibble
    
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
    
    # Now call dittoScatterPlot
    
    # Make dataframes and plot
    p.df <- dittoScatterPlot(
        object, xdat, ydat, var, shape.by, split.by,
        extra.vars, cells.use,
        show.others, size, opacity, color.panel, colors,
        split.nrow, split.ncol, NA, NA, NA, NA, NA, NA,
        assay, slot, adjustment, assay, slot, adjustment, swap.rownames,
        shape.panel, rename.var.groups, rename.shape.groups,
        min.color, max.color, min, max, order,
        xlab, ylab, main, sub, theme,
        do.hover, hover.data, hover.assay, hover.slot, hover.adjustment,
        do.contour, contour.color, contour.linetype,
        add.trajectory.lineages, add.trajectory.curves = NULL,
        trajectory.cluster.meta, trajectory.arrow.size,
        do.letter, do.ellipse, do.label, labels.size, labels.highlight,
        labels.repel, labels.split.by,
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
