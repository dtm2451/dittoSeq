################################# dittoPlot ####################################
#'
#' Plots continuous data for cutomizable cells'/samples' groupings on a y-axis
#' @import ggplot2
#'
#' @param var Single String representing a metadata or gene OR a numeric vector with length equal to the total number of cells/samples in the dataset.
#' This is the data that will be displayed.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, or the name of the object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by String representing the name of a "metadata" to use for separating the cells/samples into discrete groups. REQUIRED.
#' @param color.by String representing the name of a "metadata" to use for settin color.
#' Affects boxplot, vlnplot, and ridgeplot fills.  Default's to \code{group.by} so this input can be skipped if both are the same.
#' @param shape Integer representing a ggplot shape, OR String representing the name of a "metadata" to use for setting the shape of the jitter points.  Default = 16, dots.
#' @param cells.use String vector of cells'/samples' names which should be included OR or a logical vector that is the same length as the number of cells in the object which sets which cells to include (a.k.a. \code{USE} in \code{colnames(object)[USE]}).
#' @param plots String vectore which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot". See details section for more info.
#' @param data.type String, for when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param do.hover Logical. Default = \code{FALSE}.  If set to \code{TRUE} (and if there is a jitter plotted - the data it will work with) : object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'data.hover' argument is used to determine what data to use.
#' @param data.hover String vector, a list of variable names, c("meta1","gene1","meta2","gene2") which determines what data to show upon hover when do.hover is set to \code{TRUE}.
#' @param color.panel String vector which the set of colors to draw from.  MYcolors by default.
#' @param colors Integer vector, the indexes / order, of colors from color.panel to actually use
#' @param main String, sets the plot title. Default = "make" and if left as make, a title will be automatically generated.  To remove, set to \code{NULL}.
#' @param sub String, sets the plot subtitle
#' @param theme ggplot theme. Allows setting of a theme. Default = theme_classic when nothing is provided.
#' @param xlab String, sets the grouping-axis label (=x-axis for box and violin plots, y-axis for ridgeplots).
#' Default is \code{group.by} so it defaults to the name of the grouping information.
#' Set to \code{NULL} to remove.
#' @param ylab String, sets the continuous-axis label (=y-axis for box and violin plots, x-axis for ridgeplots).
#' Defaults to "\code{var}" or "\code{var} expression" if var is a gene.
#' @param y.breaks Numeric vector a set of breaks that should be used as major gridlines. c(break1,break2,break3,etc.) NOTE: The low and highs of this variable will override `min` and `max`.
#' @param min,max Scalars which set a custom minimum / maximum y-value to show.  Default = set based on the limits of the data in var.
#' @param x.labels String vector, c("label1","label2","label3",...) which overrides the names of the samples/groups.  NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param x.reorder Integer vector. A sequence of numbers, from 1 to the number of groupings, for rearranging their order.
#' Method: Make a first plot without this input.
#' Then treating the leftmost grouping as index 1, and the rightmost as index n.
#' Values of x.reorder should be these indices, but in the order that you would like them rearranged to be.
#' @param rotate.labels Logical. whether the labels should be rotated.  Default = \code{FALSE} = vertical labels.
#' @param add.line numeric value(s) where a dashed horizontal line should go
#' @param line.linetype String which sets the type of line.  Any ggplot linetype should work.  Defaults to "dashed"
#' @param line.color String that sets the color(s) of the horizontal line(s)
#' @param jitter.size Scalare which sets the size of the jitter shapes.
#' @param jitter.width Scalar that sets the width/spread of the jitter in the x direction. Ignored in ridgeplots.
#' @param jitter.color String which sets the color of the jitter shapes
#' @param jitter.shapes Integer vector representing which shapes to use. Default is the first in the list, 16, which corresponds to dots.
#' When jitter shape is determined by a variable, more of this input is used.
#' @param jitter.shape.legend.size Scalar which changes the size of the shape key in the legend.
#' If set to \code{NA}, \code{jitter.size} is used.
#' @param jitter.shape.legend.show Logical which sets whether the shapes lagend will be shown is \code{shape} when jitter shape is determined by a variable.
#' @param boxplot.width Scalar which sets the width/spread of the boxplot in the x direction
#' @param boxplot.color String which sets the color of the lines of the boxplot
#' @param boxplot.show.outliers Logical, whether outliers should by including in the boxplot.
#' Default is \code{FALSE} when there is a jitter plotted, \code{TRUE} if there is no jitter.
#' @param boxplot.fill Logical, whether the boxplot should be filled in or not.
#' @param vlnplot.lineweight Scalar which sets the thickness of the line that outlines the violin plots.
#' @param vlnplot.width Scalar which sets the width/spread of the jitter in the x direction
#' @param vlnplot.scaling String which sets how the widths of the of violin plots are set in relation to eachother.
#' Options are "area", "count", and "width". If the deafult is not right for your data, I recommend trying "width".
#' For a detailed explanation of each, see \code{\link{geom_violin}}.
#' @param ridgeplot.lineweight Scalar which sets the thickness of the ridgeplot outline.
#' @param ridgeplot.scale Scalar which sets the distance/overlap between ridgeplots.
#' A value of 1 means the tallest density curve just touches the baseline of the next higher one.
#' Higher numbers lead to greater overlap.  Default = 1.25
#' @param legend.show Logical. Whether the legend should be displayed. Default = \code{TRUE}.
#' @param legend.title String or \code{NULL}, sets the title for the main legend. It is set to \code{NULL} (off) by default.
#' @param data.out Logical that sets whether just the plot should be output, or a list containing the plot (\code{p}) and data (\code{data}).  Note: plotly output is turned off in this setting, but hover.data is still calculated.
#' @param ... arguments passed to dittoPlot by dittoRidgePlot.  Options are all the ones above.
#' @return a ggplot or plotly where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#' Alternatively, will return the data that would go into such a plot as well with \code{data.out=TRUE}
#' @details
#' If \code{do.hover=TRUE}, cell/sample information, determined by the \code{hover.data} input, is retrieved and displayed
#' upon hovering the cursor over a jitter point.
#'
#' If \code{data.out=TRUE}, a list containing the plot (\code{p}), a data.table containing the underlying data for target cells (\code{data})
#'
#' If \code{add.line}, is provided a value, a vertical (or horizontal, depending on plot type) line will be added at the value provided.
#' It's type and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#'
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Each plot type has specific associated options which are controlled by variables that start with their associated string, ex: \code{jitter.size}.
#'
#' Inclusion of \code{"ridgeplot"} overrides boxplot and violin plot and changes the plot to be horizontal.
#'
#' Use of dittoRidgePlot essntially just changes the default for the \code{plots} input to be \code{c("ridgeplot", "jitter")}.
#'
#' Similarly, use of dittoBoxPlot essntially just changes the default for the \code{plots} input to be \code{c("jitter", "boxplot")}.
#'
#' @examples
#' library(Seurat)
#' pbmc <- pbmc_small
#' dittoPlot("CD14", object = "pbmc", group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' dittoPlot("CD14", group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#'
#' # We can adjust the types of plots displayed with the plots input:
#' dittoPlot("CD14", group.by = "RNA_snn_res.1",
#'   plots = c("vlnplot", "boxplot", "jitter"),
#'   boxplot.fill = FALSE)
#'
#' # Quickly make a Ridgeplot
#' dittoRidgePlot("CD14", group.by = "RNA_snn_res.1")
#'
#' # Quickly make a Boxplot
#' dittoBoxPlot("CD14", group.by = "RNA_snn_res.1")
#' @export

dittoPlot <- function(var, object = DEFAULT, group.by, color.by = group.by,
                      shape = 16, cells.use = NULL, plots = c("jitter","vlnplot"),
                      data.type = "normalized", do.hover = FALSE, data.hover = var,
                      color.panel = MYcolors, colors = seq_along(color.panel),
                      theme = theme_classic(), main = "make", sub = NULL,
                      ylab = "make", y.breaks = NULL, min = NULL, max = NULL,
                      xlab = group.by, x.labels = NULL, rotate.labels = NA,
                      x.reorder = seq_along(meta.levels(group.by, object)),
                      jitter.size=1, jitter.width=0.2, jitter.color = "black",
                      jitter.shapes=c(16,15,17,23,25,8), jitter.shape.legend.size = NA,
                      jitter.shape.legend.show = TRUE,
                      boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = NA, boxplot.fill =TRUE,
                      vlnplot.lineweight = 1, vlnplot.width = 1, vlnplot.scaling = "area",
                      ridgeplot.lineweight = 1, ridgeplot.scale = 1.25,
                      add.line = NULL, line.linetype = "dashed", line.color = "black",
                      legend.show = TRUE, legend.title = NULL, data.out = FALSE){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){ object <- deparse(substitute(object)) }
  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  #Parse Title Defaults
  #ylab
  if (!(is.null(ylab)) && length(var)==1 && is.character(var)) {
      if (is.gene(var, object) && ylab == "make") { ylab <- paste0(var," expression") }
      if (ylab == "make" | ylab=="var") { ylab <- var }
  } else if (ylab == "make") { ylab <- NULL }
  #main
  if (!is.null(main) && main == "make") {
      if (length(var)==1) { main <- var } else { main <- NULL }
  }

  #Grab the data
  if (is.meta(shape, object)) { extra.vars = shape } else { extra.vars = NULL }
  Target_data <- dittoSingleAxisDataGather(main.var = var, object = object, group.by = group.by, color.by = color.by,
                                    extra.vars = extra.vars, cells.use = cells.use, data.type = data.type,
                                    do.hover = do.hover, data.hover = data.hover)$Target_data
  #Rename and/or reorder x groupings (steps 1 and 2)
  rename.args <- list(x = as.character(Target_data$grouping))
  if (!(is.null(x.reorder))){rename.args$levels <- levels(factor(rename.args$x))[x.reorder]}
  if (!(is.null(x.labels))){rename.args$labels <- x.labels}
  Target_data$grouping <- do.call(factor, args = rename.args)

  #####Start making the plot
  p <- ggplot(Target_data, aes(fill=color)) +
    theme + scale_fill_manual(name = legend.title, values=color.panel[colors]) +
    ggtitle(main, sub)
  #Add data
  if(!("ridgeplot" %in% plots)) {
    p <- dittoYPlotMaker(
      p, Target_data, plots, xlab, ylab, shape, jitter.size, jitter.width,
      jitter.color,jitter.shapes, jitter.shape.legend.size, jitter.shape.legend.show,
      boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
      vlnplot.lineweight, vlnplot.width, vlnplot.scaling, add.line,
      line.linetype, line.color, rotate.labels, do.hover, y.breaks, min, max)
  } else {
    p <- dittoXPlotMaker(
      p, Target_data, plots, xlab, ylab, jitter.size, jitter.color,
      jitter.shape.legend.size, jitter.shape.legend.show,
      ridgeplot.lineweight, ridgeplot.scale, add.line, line.linetype,
      line.color, rotate.labels, do.hover, color.panel, colors, y.breaks, min, max)}
  #Remove legend, if warrented
  if (!legend.show) { p <- remove_legend(p) }
  #DONE. Return the plot
  if(data.out) {return(list(p = p, data = Target_data))}
  else { if(do.hover & ("jitter" %in% plots)){
    return(plotly::ggplotly(p, tooltip = "text")) }
    else { return(p) }
  }
}

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings horizontally in a desnsity representation
#' @export
dittoRidgePlot <- function(..., plots = c("ridgeplot")){ dittoPlot(..., plots = plots) }

#' @describeIn dittoPlot Plots continuous data for cutomizable cells'/samples' groupings in boxplot form
#' @export
dittoBoxPlot <- function(..., plots = c("boxplot","jitter")){ dittoPlot(..., plots = plots) }

dittoYPlotMaker <- function(
  p, Target_data, plots, xlab, ylab, shape, jitter.size, jitter.width,
  jitter.color,jitter.shapes, jitter.shape.legend.size, jitter.shape.legend.show,
  boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
  vlnplot.lineweight, vlnplot.width, vlnplot.scaling, add.line,
  line.linetype, line.color, rotate.labels, do.hover, y.breaks, min, max) {
  # This function takes in a partial dittoPlot ggplot object without any data overlay,
  # and parses adding the main data visualizations.
  # It adds plots based on what is requested in plots, ordered by their order.

  # Now that we know the plot's direction, set y-axis limits
  if (!is.null(y.breaks)) {
    p <- p + scale_y_continuous(breaks = y.breaks) + coord_cartesian(ylim=c(min(y.breaks),max(y.breaks)))
  } else {
    if (is.null(min)){min <- min(Target_data$var.data)}
    if (is.null(max)){max <- max(Target_data$var.data)}
    p <- p + coord_cartesian(ylim=c(min,max))
  }

  # Add labels and, if requested, lines
  p <- p + aes(x = grouping, y = var.data) + xlab(xlab) + ylab(ylab)
  if (is.na(rotate.labels)) {rotate.labels <- TRUE}
  if (rotate.labels) {p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1, size=12))}
  if (!is.null(add.line)) {p <- p + geom_hline(yintercept=add.line, linetype= line.linetype, color = line.color)}

  # Add Plots
  for (i in seq_along(plots)){
    if (plots[i] == "boxplot") {
      boxplot.args <- list(width=boxplot.width, color = boxplot.color, alpha = ifelse(boxplot.fill, 1, 0))
      if (is.na(boxplot.show.outliers)){ boxplot.show.outliers <- ifelse("jitter" %in% plots, FALSE, TRUE) }
      if (!boxplot.show.outliers){ boxplot.args$outlier.shape = NA }
      p <- p + do.call(geom_boxplot, boxplot.args)
    }
    if (plots[i] == "jitter") {
      jitter.args <- list(size=jitter.size, width=jitter.width, height = 0, color = jitter.color)
      #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
      if (is.character(shape)){
        #Make jitter with shapes
        jitter.args$mapping <- if (do.hover) {aes(shape = shape, text = hover.string)} else {aes(shape = shape)}
        p <- p + do.call(geom_jitter, jitter.args) +
          scale_shape_manual(
            values = jitter.shapes[seq_along(levels(as.factor(Target_data[,shape==names(Target_data)])))],
            labels = levels(as.factor(as.character(Target_data[,shape==names(Target_data)]))))
        if (!is.na(jitter.shape.legend.size)){
          p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
        }
        if (jitter.shape.legend.show==FALSE){ p <- p + guides(shape = "none") }
      } else {
        jitter.args$mapping <- if (do.hover) {aes(text = hover.string)}
        p <- p + do.call(geom_jitter, jitter.args)
      }
    }
    if (plots[i] == "vlnplot") {
      p <- p + geom_violin(size = vlnplot.lineweight, width = vlnplot.width, scale = vlnplot.scaling)
    }
  }
  p
}

#' @importFrom ggridges geom_density_ridges2
dittoXPlotMaker <- function(p, Target_data, plots, xlab, ylab,
                            jitter.size=1, jitter.color = "black",
                            jitter.shape.legend.size, jitter.shape.legend.show,
                            ridgeplot.lineweight = 1, ridgeplot.scale = 1.25,
                            add.line=NULL, line.linetype = "dashed", line.color = "black",
                            rotate.labels = FALSE, do.hover, color.panel, colors, y.breaks, min, max) {
  #This function takes in a partial dittoPlot ggplot object without any data overlay, and parses adding the main data visualizations.
  # It adds plots based on what is requested in plots, *ordered by their order*

  # Now that we know the plot's direction, set continuous-axis limits
  # if (!is.null(y.breaks)) {
  #   p <- p + scale_x_continuous(breaks = y.breaks, expand = expand_scale(mult=c(0, 0))) +
  #     coord_cartesian(xlim=c(min(y.breaks),max(y.breaks)))
  # } else {
  #   if (is.null(min)){
  #     min <- min(Target_data$var.data)
  #     exp.min <- 2} else {exp.min = 1}
  #   if (is.null(max)){
  #     max <- max(Target_data$var.data)
  #     exp.max <- 2} else {exp.max = 1}
  #   p <- p + coord_cartesian(xlim=c(min,max)) +
  #   scale_x_continuous(expand = expand_scale(mult=c(exp.min, exp.max)))
  # }

  # Add labels and, if requested, lines
  p <- p + aes(x = var.data, y = grouping) + xlab(ylab) + ylab(xlab) +
  scale_y_discrete(expand = expand_scale(mult = c(0.02, 0.05)))

  if (is.na(rotate.labels)) {rotate.labels <- FALSE}
  if (rotate.labels) {p <- p + theme(axis.text.y= element_text(angle=45, hjust = 1, vjust = 1, size=12))}
  if (!is.null(add.line)) {p <- p + geom_vline(xintercept=add.line, linetype= line.linetype, color = line.color)}
  p <- p + ggridges::geom_density_ridges2(
    size = ridgeplot.lineweight, scale = ridgeplot.scale,
    jittered_points = "jitter" %in% plots, point_size = jitter.size, point_color = jitter.color) +
    scale_color_manual(values=color.panel[colors])
  if (!is.na(jitter.shape.legend.size)){
    p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
  }
  if (jitter.shape.legend.show==FALSE){ p <- p + guides(shape = "none") }
  p
}

#' Gathers data for single-axis plotting functions
#' @import ggplot2
#'
#' @param main.var Character representing the name of a metadata, gene, or "ident" to be plotted. OR character, numeric, or factor vector of length equal to the total # cells in the dataset.
#' @param object                 the Seurat or RNAseq Object = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by               "metadata" to use for separating values. REQUIRED.
#' @param color.by               "metadata" to use for coloring. Affects boxplot and vlnplot fills. REQUIRED when using either.
#' @param extra.vars               "metadata" to use for setting the shape of jitter.  Default = just dots. Ignored if not a "quoted" metadata or "ident"
#' @param cells.use              Cells to include: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. USE in object@cell.names[USE])
#' @param data.type              For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param do.hover               TRUE/FALSE. Default = FALSE.  If set to true (and if there is a jitter plotted - the data it will work with) : object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'data.hover' argument is used to determine what data to use.
#' @param data.hover             list of variable names, c("meta1","gene1","meta2","gene2"). determines what data to show on hover when do.hover is set to TRUE.
#' @return Generates Target_data and Others_data data.frames for use by single-axis DittoSeq plotters.

dittoSingleAxisDataGather <- function(main.var, object = DEFAULT, group.by = "Sample", color.by = group.by,
                                      extra.vars = NULL, cells.use = NULL, data.type = "normalized",
                                      do.hover = FALSE, data.hover = c(main.var, extra.vars)){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){ object <- deparse(substitute(object)) }
  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)
  ###Make dataframe for storing the plotting data:
  full_data <- data.frame(var.data = var_OR_get_meta_or_gene(main.var, object, data.type),
                          grouping = meta(group.by, object),
                          color = meta(color.by, object),
                          row.names = all.cells)
  names <- names(full_data)
  #Add Extra data
  if(length(extra.vars)>0){
    for (i in seq_along(extra.vars)){
      full_data <- cbind(full_data, var_OR_get_meta_or_gene(extra.vars, object, data.type))
    }
    names <- c(names, extra.vars)
  }
  #Add hover strings
  if (do.hover) {
    full_data$hover.string <- make_hover_strings(data.hover, object, data.type)
    names <- c(names, data.hover)
  }
  colnames(full_data) <- names

  return(list(Target_data = full_data[all.cells %in% cells.use,],
              Others_data = full_data[!(all.cells %in% cells.use),]))
}

