#' Plot spatial transcriptomics data by spatial coordinates
#' @inheritParams dittoDimPlot
#' @param split.by 1 or 2 strings naming discrete metadata to use for splitting the cells/samples into multiple plots with ggplot faceting.
#' *Defaults here to \code{samples.meta}, but otherwise works the ssame as in other functions:
#'
#' When 2 metadata are named, c(row,col), the first is used as rows and the second is used for columns of the resulting grid.
#'
#' When 1 metadata is named, shape control can be achieved with \code{split.nrow} and \code{split.ncol}
#'
#' @param samples.use String vector (or NULL meaning 'use all') for subsetting to certain tissue slices.
#' It's elements should be the options of \code{samples.meta} values that should be retained. 
#' @param cells.use An additional method for subsetting cells to show, the retained cells/locations/spots will be in intersect of those indicated by this input & \code{samples.use}.
#' Otherwise, this input works exactly the same here as in other dittoSeq functions:
#' 
#' String vector of cells'/samples' names OR an integer vector specifying the indices of cells/samples which should be included.
#' 
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' 
#' @param samples.meta String naming the metadata that contains tissue slice (+ image) linkage information.
#' @param images.id Single string, a value of \code{imgData(object)$image_id} when \code{object} is a \code{\link[SpatialExperiment]{VisiumExperiment}}, for selecting a set of tissue images to use.
#' @param image.paths Alternative to \code{images.id} which takes precedence when given.
#' A single string, or string vector, specifying file path(s) of tissue image(s).
#' If a particular tissue slice should not have an image displayed, use the value NA.
#' Order is assumed to be the same as that of \code{samples.use}.
#' @return a ggplot scatterplot where colored dots and/or shapes represent data from the given spatial coordinates.
#'
#' Alternatively, if \code{data.out=TRUE}, a list containing three slots is output: the plot (named 'p'), a data.table containing the underlying data for target cells (named 'Target_data'), and a data.table containing the underlying data for non-target cells (named 'Others_data').
#'
#' Alternatively, if \code{do.hover} is set to \code{TRUE} and no images are to be shown (incompatibility), the plot is converted from ggplot to plotly &
#' cell/sample information, determined by the \code{hover.data} input, is retrieved, added to the dataframe, and displayed upon hovering the cursor over the plot.
#'
#' @details X/Y are extracted... Remember to update notes on defaults for these inputs: show.gridlines
#' @export
#' @author Daniel Bunis


dittoSpatial <- function(
    object,
    var = NULL,
    samples.meta = "sample_id",
    samples.use = NULL,
    cells.use = NULL,
    split.by = samples.meta,
    images.id = "lowres",
    image.paths = NULL,
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
    show.grid.lines = FALSE,
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
    # add.trajectory.curves = NULL,
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
        xdat <- coords$array_row*SpatialExperiment::scaleFactors(object)
        ydat <- coords$array_col*SpatialExperiment::scaleFactors(object)
    } else {
        xdat <- coords$x
        ydat <- coords$y
    }

    ### Combine all cell subsetting inputs into cells.use
    # Interpret samples.meta & subset by samples.use, intersect with cells.use
    cells.use <- .which_cells(cells.use, object)
    if (!is.null(samples.use)) {
        cells.in.samples <- 
            .which_cells(meta(samples.meta, object) %in% samples.use, object)
        cells.use <- intersect(cells.use, cells.in.samples)
    } else {
        samples.use <- metaLevels(samples.use, object, cells.use)
    }
    
    ### Interpret images.id & images.paths
    grobs <- NULL
    
    # If paths given, use those. (Assume linked to samples.use order if more than one.)
    if (!is.null(image.paths)) {
        grobs <- lapply(
            image.paths,
            function(x) {
                if (is.na(x)) {
                    NA
                } else {
                    img <- magick::image_read(x)
                    grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
                }
            })
    } else if (is_visium && !is.null(images.id)) {
        # Determine relevant images per sample if given as an id (default).
        img_data <- SpatialExperiment::imgData(object)
        img_data[img_data$sample_id %in% samples.use,]
        
        grobs <- lapply(
            seq_len(nrow(img_data))[img_data$image_id == images.id],
            function(img_ind) {
                if (!is.null(SpatialExperiment::imgGrob(img_data$data[img_ind]))) {
                    SpatialExperiment::imgGrob(img_data$data[img_ind])
                } else if (SpatialExperiment::imgPath(img_data$data[img_ind])) {
                    img <- magick::image_read(SpatialExperiment::imgPath(img_data$data[img_ind]))
                    grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
                } else {
                    message("imgData for ", img_data$sample_id[img_ind], "not found or not supported.")
                    NA
                }
        })
    }
    
    # Turn image grobs into tibble for easy ggplot compatibility.
    if (!is.null(grobs)) {
        images_tibble <- dplyr::tibble(
            height = lapply(grobs, function(x) { nrow(x$raster) }),
            width = lapply(grobs, function(x) { ncol(x$raster) }),
            grob = grobs
        )
        # Add split.by unless single image
        if (length(grobs)>1) {
            images_tibble[[samples.meta]] <- samples.use
        }
        # Remove NA grobs
        images_tibble <- images_tibble[!is.na(images_tibble$grob),]
    }
    
    if (nrow(images_tibble) > 1) {
        images_tibble
    }
    
    
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
    
    # Make dataframes and plot with dittoScatterPlot
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
