#' Show RNAseq data overlayed on a scatter plot
#' @import ggplot2
#'
#' @param x.var Variable for setting x-axis position of cells/samples.  Can be the name of a gene, meta-data, or "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param y.var Variable for setting y-axis position of cells/samples.  Can be the name of a gene, meta-data, or "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param overlay.color.var Variable for setting the color of cells/samples in the plot.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param overlay.shape.var Variable for setting the shape of cells/samples in the plot.  Note: must be discrete.  Can be the name of a gene or meta-data.  Alternatively, can be "ident" for clusters of a Seurat object.  Alternatively, can be a numeric of length equal to the total number of cells/samples in object.
#' @param object the Seurat, SingleCellExperiment, or RNAseq object to work on.
#' @param cells.use cells to show: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. *THIS*: object@cell.names[*THIS*])
#' @param show.others TRUE/FALSE. TRUE by default, whether other cells should be shown in the background
#' @param color.panel a list of colors to be used for when plotting a discrete var.
#' @param colors indexes / order of colors from color.panel to use. USAGE= changing the order of how colors are linked to specific groups
#' @param data.type.x For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param data.type.y For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param data.type.color For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized"
#' @param do.hover TRUE/FALSE. Default = FALSE.  If set to true: object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'data.hover' argument is used to determine what data to use.  NOTE: incompatible with lettering (due to a ggplotly incompatibility). Setting do.hover to TRUE will override a do.letter=TRUE or NA.
#' @param data.hover list of variable names, c("meta1","gene1","meta2","gene2"). determines what data to show on hover when do.hover is set to TRUE.
#' @param data.type.hover For when plotting expression data, sets the data-type slot that will be obtained. See \link[dittoSeq]{gene}. DEFAULT = "normalized" For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param shape Number for setting the shape
#' @param shapes If multiple shapes are needed, this sets the sahpes to pull from.  Default is a list of 6, \code{c(16,15,17,23,25,8)}.  There are more, but not many of the default ggplot options are great.  Unfortunately, shape is also harder to see when points are on top of each other.  For these reasons, even as a color blind person myself, I recommend use of colors for variables with 7+ options.
#' @param size Number. Size of data points.  Default = 1.
#' @param opacity Number between 0 and 1. Great for when you have MANY overlapping points, this sets how see-through the points should be; 1 = not at all; 0 = invisible. Default = 1.
#' @param rename.color.groups Character vector containing new names for the identities of the color overlay.
#' @param rename.shape.groups Character vector containing new names for the identities of the shape overlay.
#' @param legend.show TRUE/FALSE. Whether the legend should be displayed. Default = TRUE.
#' @param legend.color.title For adding a title to the colors legend.  Set to \code{NULL} to turn off.
#' @param legend.color.size Override for legend size of the color variable.
#' @param legend.shape.title For adding a title for the colors legend.  Set to \code{NULL} to turn off
#' @param legend.shape.size Override for legend size of the shape variable.
#' @param min.color color for lowest values shown of the color overlay.  Default is \code{dittoColors()[4]}, a yellow.
#' @param max.color color for highest values shown of the color overlay.  Default is \code{dittoColors()[5]}, a dark blue.
#' @param min set the value associated with the minimum color.  All points with a lower value than this will get the same min.color.
#' @param max set the value associated with the maximum color.  All points with a higher value than this will get the same max.color.  Note: if your legend is not plotting, it may be because min > max.
#' @param breaks Numeric vector. Sets the discrete values to show in the color-scale legend for continuous data.
#' @param breaks.labels String vector with same length as \code{breaks}. Renames the values displayed next to the color-scale.
#' @param main plot title.  Default = \code{NULL}
#' @param sub plot subtitle.  Default = \code{NULL}
#' @param xlab label for y axes.  Default = \code{x.var}. To remove, set to NULL.
#' @param ylab label for y axes.  Default = \code{x.var}. To remove, set to NULL.
#' @param theme Allows setting of a theme. Default = theme_bw when nothing is provided.
#' @param data.out Whether just the plot should be output, or a list with the plot and Target_data and Others_data dataframes.  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @return Makes a plot where colored dots and/or shapes representing individual cells/samples are overlayed onto a scatterplot where x and y can be gene expression (or any numeric metadata) of those cells/samples.
#' @export
#' @examples
#'
#' library(Seurat)
#' pbmc <- pbmc_small
#' dittoScatterPlot(
#'     x.var = "nCount_RNA", y.var = "nFeature_RNA",
#'     object = "pbmc", overlay.color.var = "RNA_snn_res.1")
#'
#' # Note: scatterplots like this can be very useful for dataset QC, epecially
#' #   with percentage of reads coming from genes as the color overlay.
dittoScatterPlot <- function(x.var, y.var, overlay.color.var = NULL, overlay.shape.var = NULL,
                             object = DEFAULT, cells.use = NULL, show.others = FALSE,
                             color.panel = dittoColors(), colors = seq_along(color.panel),
                             data.type.x = "normalized", data.type.y = "normalized",
                             data.type.color = "normalized", do.hover = FALSE, data.hover = NULL, data.type.hover = "normalized",
                             shape = 16, shapes=c(16,15,17,23,25,8), size = 1, opacity = 1,
                             rename.color.groups = NA, rename.shape.groups = NA,
                             legend.show = TRUE,
                             legend.color.title = overlay.color.var, legend.color.size = 5,
                             legend.shape.title = overlay.shape.var, legend.shape.size = 5,
                             min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
                             breaks = waiver(), breaks.labels = waiver(),
                             xlab = x.var, ylab = y.var, main = NULL, sub = NULL, theme = theme_bw(),
                             data.out = FALSE){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){
    object <- deparse(substitute(object))
  }
  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- .which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- .all_cells(object)
  #Grab data
  dat <- data.frame(X = .var_OR_get_meta_or_gene(x.var, object, data.type.x),
                    Y = .var_OR_get_meta_or_gene(y.var, object, data.type.y),
                    row.names = all.cells)
  if(!(is.null(overlay.color.var))){
    dat$color <- .var_OR_get_meta_or_gene(overlay.color.var, object, data.type.color)
    do.color <- TRUE
    is.color.numeric <- is.numeric(dat$color)
  } else {do.color = FALSE}
  if(!(is.null(overlay.shape.var))){
    dat$shape <- .var_OR_get_meta_or_gene(overlay.shape.var, object, data.type.color)
  }
  #Grab data for hovering
  if (do.hover) {
    hover.string <- .make_hover_strings_from_vars(data.hover, object, data.type.hover)
  } else {hover.string <- NA}
  #Split data to target vs other
  Target_data <- dat[cells.use,]
  Others_data <- dat[!(all.cells %in% cells.use),]
  ###Start building the plot###
  p <- ggplot() + ylab(ylab) + xlab(xlab) + ggtitle(main,sub) + theme
  if(!(is.null(overlay.shape.var))){ p <- p +
    scale_shape_manual(values = shapes[seq_along(levels(as.factor(Target_data$shape)))],
                       label = if (!(is.na(rename.shape.groups[1]))){rename.shape.groups}
                               else {levels(as.factor(Target_data$shape))},
                       name = legend.shape.title) +
    guides(shape = guide_legend(override.aes = list(size=legend.shape.size)))
  }
  if(do.color){
    if (is.color.numeric){ p <- p +
      scale_colour_gradient(low= min.color, high = max.color,
                            limits = c(ifelse(is.null(min), min(Target_data$color), min),
                                       ifelse(is.null(max), max(Target_data$color), max)),
                            breaks = breaks, labels = breaks.labels,
                            name = legend.color.title)
    } else {
      args.colour <- list(
          name = legend.color.title,
          values = color.panel[colors])
      if (!(is.na(rename.color.groups[1]))) {
          args.colour$label <- rename.color.groups
      }
      p <- p +
          do.call(scale_colour_manual, args.colour) +
          guides(color = guide_legend(override.aes = list(size=legend.color.size)))
    }
  }
  ###Add the data###
  #Make gray dots on the bottom layer if show.others = TRUE and cells.use is a subset of all the cells / samples.
  if (show.others & dim(Others_data)[1]>1) {
    p <- p + geom_point(data=Others_data,
                        if(do.hover){aes(x = X, y = Y, text = hover.string)}else{aes(x = X, y = Y)},
                        size=size, color = "gray90")
  }
  #Overlay the target data on top
  # If 'shape' input was the name of a meta.data, aka type=character, treat shape as an aesthetic for performing grouping.
  # Otherwise it is a number and belongs outside of aes.
  if(!(is.null(overlay.color.var)) & !(is.null(overlay.shape.var))){
    p <- p + geom_point(data=Target_data,
                        if(do.hover){aes(x = X, y = Y, colour = color, shape=shape, text = hover.string)}
                        else{aes(x = X, y = Y, colour = color, shape=shape)},
                        size=size, alpha = opacity)
  }
  if(!(is.null(overlay.color.var)) & is.null(overlay.shape.var)){
    p <- p + geom_point(data=Target_data,
                        if(do.hover){aes(x = X, y = Y, colour = color, text = hover.string)}
                        else{aes(x = X, y = Y, colour = color)},
                        shape= shape, size=size, alpha = opacity)
  }
  if(is.null(overlay.color.var) & !(is.null(overlay.shape.var))){
    p <- p + geom_point(data=Target_data,
                        if(do.hover){aes(x = X, y = Y, shape=shape, text = hover.string)}
                        else{aes(x = X, y = Y, shape=shape)},
                        colour= "black", size=size, alpha = opacity)
  }
  if(is.null(overlay.color.var) & is.null(overlay.shape.var)){
    p <- p + geom_point(data=Target_data,
                        if(do.hover){aes(x = X, y = Y, text = hover.string)}
                        else{aes(x = X, y = Y)},
                        colour= "black", shape = shape, size=size, alpha = opacity)
  }
  #Remove legend, if warrented
  if (!legend.show) { p <- remove_legend(p) }
  ### RETURN the PLOT ###
  if(data.out){return(list(plot = p,
                           Target_data = Target_data,
                           Others_data = Others_data))}
  else{
    if(do.hover){ return(plotly::ggplotly(p, tooltip = "text")) }
    else { return(p) }
  }
}
