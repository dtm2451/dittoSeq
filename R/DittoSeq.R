################# DBDimPlot ####################

#' Show data overlayed on a tsne or pca or other reduction-type
#'
#' @param var Target Variable = either values or a metadata (in "quotes"), gene (in "quotes"), or "ident"
#' @param object the Seurat or RNAseq object to work on
#' @param reduction.use "pca", "tsne", "ica", etc.  Default = tsne for Seurat objects, and pca for RNAseq objects
#' @param dim.1 The component number to use on the x-axis.  Default = 1
#' @param dim.2 The component number to use on the y-axis.  Default = 2
#' @param theme Allows setting of a theme. Default = theme_bw when nothing is provided.
#' @param size Number. Size of data points.  Default = 1.
#' @param shape Number for setting shape OR name of metadata to use for setting shape
#' @param shapes the shapes to use.  Default is a list of 6.  There are more, but not many of the default ggplot options are great.  I recommend using colors for any variable with 7+ options.
#' @param legend.size The size to increase the plotting of colors/letters legend shapes to (for discrete variable plotting)
#' @param legend.title For adding a title for the colors/letters legend.  It is set to NULL (off) by default.
#' @param shape.legend.size The size to increase the plotting of shapes legend shapes to.
#' @param shape.legend.title For adding a title for the shapes legend is a meta.data was given to 'shape' and multiple shapes were therefore used.  It is set to NULL (off) by default.
#' @param data.type For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param main plot title
#' @param sub plot subtitle
#' @param xlab label for y axes.  Default labels are generated if you do not give this a specific value.  To remove the labeling, set to NULL.
#' @param ylab label for y axes.  Default labels are generated if you do not give this a specific value.  To remove the labeling, set to NULL.
#' @param auto.title TRUE/FALSE = whether a default 'main' should be generated.
#' @param cells.use cells to show: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. *THIS*: object@cell.names[*THIS*])
#' @param show.others TRUE/FALSE. TRUE by default, whether other cells should be shown in the background
#' @param ellipse TRUE/FALSE. Whether the groups should be surrounded by an ellipse.
#' @param do.label  TRUE/FALSE. Whether to add text labels at the center (median) of clusters for grouping vars
#' @param label.size Size of the the labels text
#' @param highlight.labels TRUE/FALSE. Whether the labels should have a box behind them
#' @param rename.groups new names for the identities of var.  Change to NULL to remove labeling altogether.
#' @param min.color color for lowest values of var/min
#' @param max.color color for highest values of var/max
#' @param min set the value associated with the minimum color.  All points with a lower value than this will get the same min.color.
#' @param max set the value associated with the maximum color.  All points with a higher value than this will get the same max.color.  Note: if your legend is not plotting, it's likely because min > max.
#' @param color.panel a list of colors to be used for when plotting a discrete var.
#' @param colors indexes / order of colors from color.panel to use. USAGE= changing the order of how colors are linked to specific groups
#' @param do.letter TRUE/FALSE/NA. Whether letters should be added on top of the colored dots. For colorblindness compatibility.  NA by default, and if left that way, will be set to either TRUE or FALSE depending on the number of groups if a discrete var is given. For when there are lots of descrete variables, it can be hard to see by just color / shape.  NOTE: Lettering is incompatible with changing dots to shapes, so this will override a setting of distinct shapes based on the 'shape' variable!
#' @param do.hover TRUE/FALSE. Default = F.  If set to true: object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'data.hover' argument is used to determine what data to use.  NOTE: incompatible with lettering (due to a ggplotly incompatibility). Setting do.hover to TRUE will override a do.letter=TRUE or NA.
#' @param data.hover list of variable names, c("meta1","gene1","meta2","gene2"). determines what data to show on hover when do.hover is set to TRUE.
#' @return Makes a plot where colored dots (or other shapes) are overlayed onto a tSNE, PCA, ICA, ..., plot of choice.  var is the argument that sets how dots will be colored, and it can refer to either continuous (ex: "CD34" = gene expression) or discrete (ex: "ident" = clustering) data.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' DBDimPlot("res.1", object = "pbmc")
#' DBDimPlot("ident", object = "pbmc", reduction.use = "pca", ellipse = TRUE, do.label = TRUE)
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBDimPlot("res.1")
#' DBDimPlot("ident", reduction.use = "pca", ellipse = TRUE, do.label = TRUE)

DBDimPlot <- function(var="ident", object = DEFAULT, reduction.use = NA, dim.1 = 1, dim.2 = 2, theme = NA,
                      size=1, shape=16, shapes=c(16,15,17,23,25,8),
                      legend.size = 5, legend.title = NULL,
                      shape.legend.size = 5, shape.legend.title = NULL,
                      data.type = "normalized",
                      main = NULL, sub = NULL, xlab="make", ylab ="make", auto.title = T,
                      cells.use = NULL, show.others=TRUE, ellipse = F,
                      do.label = F, label.size = 5, highlight.labels = T,
                      rename.groups = NA,
                      min.color = "#F0E442", max.color = "#0072B2", min = NULL, max = NULL,
                      color.panel = MYcolors, colors = 1:length(color.panel),
                      do.letter = NA, do.hover = FALSE, data.hover = var){

  #Change object to character if not already
  object <- S4_2string(object)

  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  #If reduction.use = NA (was not provided), populate it to be tsne or pca.
  if (grepl("Seurat",classof(object)) & is.na(reduction.use)) {reduction.use <- "tsne"}
  if (classof(object)=="RNAseq" & is.na(reduction.use)) {reduction.use <- "pca"}

  #Generate the x/y dimensional reduction data and axes labels.
  xdat <- extDim(reduction.use, dim.1, object)
  ydat <- extDim(reduction.use, dim.2, object)
  #If xlab/ylab left as "make", use default axis labels generated in the extDim call (ex. "tSNE_1" or "PC2").
  if (!(is.null(xlab))) {
    if (xlab=="make") {
      xlab <- xdat$name
    }
  }
  if (!(is.null(ylab))) {
    if (ylab=="make") {
      ylab <- ydat$name
    }
  }

  #Build data for populating dat, the data.frame for plotyting.
  #Determine the identity of the provided 'var' and populate Y, the variable used for coloring.
  Y <- var_OR_get_meta_or_gene(var, object, data.type)

  #Decide if letters should be added or not
  #If#1) if do.hover was set to TRUE, lettering will not work, so set to FALSE.
  if(do.hover){
    do.letter <- FALSE
  } else {
    #IF#2) if do.letter was already set, we'll just go with what the user wanted!
    if(is.na(do.letter)){
      #If#2) if the data is discrete, continue. (Otherwise, it is continuous so letters would not make sense!)
      if(!(is.numeric(Y))){
        #If#3) if the number of groups is 8 or more, letters are recommended.
        if(length(levels(as.factor(Y)))>=8){
          do.letter <- TRUE
        } else { do.letter <- FALSE }
      } else { do.letter <- FALSE }
    }
  }
  #Determine the identity of the provided 'shape'.
  #If it is a meta.data name, pull the meta.  else (it is a number) just carry it through
  if(typeof(shape)=="character"){
    if(do.letter){Shape <- shapes[1]
    } else {Shape <- meta(shape, object)}
  } else {Shape <- shape}

  #Groundwork for plotly hover data:
  #This does: if do.hover=T and data.hover has a list of genes / metas,
  # then for all cells, make a string "var1: var1-value\nvar2: var2-value..."
  hover.string <- NA
  if (do.hover) {
    hover.string <- make_hover_strings(data.hover, object, data.type)
  }

  #Populate the data.frame to be used for plotting
  full_dat <- data.frame(Y = Y,
                         dim1 = xdat$embeddings,
                         dim2 = ydat$embeddings,
                         size = size,
                         shape = Shape,
                         Group = Y,
                         hover.string = hover.string)

  #Subset to cells.use
  Others_dat <- full_dat[!(all.cells %in% cells.use),]
  Target_dat <- full_dat[cells.use,]

  ## Groundwork for adding letter labeling of dots
  # If adding letters, create a vector of what those labels should be
  if(do.letter){
    letter.labels <- c(LETTERS, letters, 0:9, "!", "@", "#", "$", "%", "^", "&", "*", "(", ")", "-",
                       "+", "_", "=", ";", "/", "|", "{", "}", "~")[1:length(levels(as.factor(Y)))]
    names(letter.labels) <- levels(Y)
    letter.colors <- c(rep("white",7), rep("black",10),rep("white",9))[1:length(levels(as.factor(Y)))]
    names(letter.colors) <- levels(Y)
  }

  ###Start building the plot###
  p <- ggplot() + ylab(ylab) + xlab(xlab)

  #Then Add more layers:

  ###Add the data###
  #Make gray dots on the bottom layer if show.others = T and cells.use is a subset of all the cells / samples.
  if (show.others & dim(Others_dat)[1]>1) {
    p <- p + geom_point(data=Others_dat,
                        if(do.hover){aes(x = dim1, y = dim2, text = hover.string)}else{aes(x = dim1, y = dim2)},
                        size=0.5, color = "gray90")
  }
  #Overlay the target data on top
  # If 'shape' input was the name of a meta.data, aka type=character, treat shape as an aesthetic for performing grouping.
  # Otherwise it is a number and belongs outside of aes.
  if (typeof(shape)=="character" & !do.letter) {
    p <- p + geom_point(data=Target_dat,
                        if(do.hover){aes(x = dim1, y = dim2, colour = Y, shape=shape, text = hover.string)
                        }else{aes(x = dim1, y = dim2, colour = Y, shape=shape)},
                        size=size) +
      # if(do.letter){geom_point(data=Target_dat, aes(x = dim1, y = dim2, shape = Y), color = letter.colors, size=size)} +
      scale_shape_manual(values = shapes[1:length(levels(as.factor(Target_dat$shape)))],
                         labels = levels(as.factor(as.character(Target_dat$shape))))
  }  else {
    p <- p + geom_point(data=Target_dat,
                        if(do.hover){aes(x = dim1, y = dim2, colour = Y, text = hover.string)}
                        else{aes(x = dim1, y = dim2, colour = Y)},
                        shape= Shape, size=size, stroke = 0)
    if(do.letter){
      p <- p +
        geom_point(data=Target_dat, aes(x = dim1, y = dim2, shape = Y), color = "black", size=size/2) +
        scale_shape_manual(name = legend.title,
                           values = letter.labels#,
                           # if(!(is.na(rename.groups[1]))){
                           #   labels = rename.groups}
        )
    }
  }

  ###Add ellipse###
  ### Draw an ellipse if ellipse = T.
  if (ellipse) { p <- p + stat_ellipse(data=Target_dat,
                                       aes(x = dim1, y = dim2, colour = Y),
                                       type = "t",
                                       linetype = 2,
                                       size = 0.5,
                                       show.legend = F
  )}

  ###Add titles###
  #If not provided, autogenerate based on the identity of var
  if (is.null(main) & auto.title==T & length(var)==1){
    #If var is a meta.data, make the title = var
    if(is.meta(var, object)){main <- var}
    #If var is a gene, still make the title = var
    if(is.gene(var, object)){main <- var}
  }
  #Add titles
  p <- p + ggtitle(main, subtitle = sub)

  ### Add Labels ###
  if (do.label) {
    #Make a text plot at the median x and y values for each cluster
    #Determine medians
    cent.1 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim1[Target_dat$Y==X]))
    cent.2 = sapply(levels(as.factor(Target_dat$Y)), function(X) median(Target_dat$dim2[Target_dat$Y==X]))
    #Add labels
    if (highlight.labels){
      #Add labels with a white background
      p <- p +
        geom_label(data = data.frame(x=cent.1, y=cent.2),
                   aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(as.factor(Target_dat$Y))}),
                   size = label.size)
    } else {
      #Add labels without a white background
      p <- p +
        geom_text(data = data.frame(x=cent.1, y=cent.2),
                  aes(x = x, y = y, label = if(!(is.na(rename.groups[1]))){rename.groups} else {levels(as.factor(Target_dat$Y))}),
                  size = label.size)
    }
  }

  ### Set the colors ###
  ### Also change the size of the dots in the legend if showing groupings ###
  #If var yielded a list of groups for plotting (should be in the form of a list of strings = character, or a factor = integer)
  if (!(is.numeric(Y))){
    #If the number of levels/groups is less than the length of the color.panel, use the color.panel set.
    if (length(levels(as.factor(as.character(Target_dat$Y))))<=length(color.panel[colors])){
      #If labels.rename was changed from NA, rename the grouping labels to rename.groups values
      if (!(is.na(rename.groups[1]))){
        p <- p+ scale_colour_manual(name = legend.title,
                                    values = color.panel[colors],
                                    labels = rename.groups)
      } else {
        #If not, just set the colors and name for the legend key
        p <- p+ scale_colour_manual(name = legend.title,
                                    values = color.panel[colors])
      }
      #Also change the legend properties.
      if (!is.null(legend.size)){
        #First, change the size of the dots in the legend, unless legend.size was set to NULL
        #Also, set the legend title.  Given input if given, or remove if still null.
        if(do.letter){
          p <- p + guides(colour = guide_legend(override.aes = list(size=legend.size)))
        } else {
          p <- p + guides(colour = guide_legend(override.aes = list(size=legend.size)),
                          shape = guide_legend(override.aes = list(size=shape.legend.size), title = shape.legend.title))
        }
      }
    }
  } else {
    #Otherwise, the data is continous, so set a gradient that goes from 'min.color' input color to 'max.color' input color.
    p <- p + scale_color_gradient(low= min.color, high = max.color, limits = c(ifelse(is.null(min), min(Target_dat$Y), min),
                                                                               ifelse(is.null(max), max(Target_dat$Y), max)),
                                  #Next, set the legend title.  Given input if given, or remove if still null.
                                  name = legend.title)
  }

  ### Set the theme ###
  #Use theme_bw if 'theme' = NA (was not provided), or use prettyplot.1 if "prettyplot" was provided, or
  # use provided theme if a full one is provided, aka = a list.
  if (is.na(theme)){
    p <- p + theme_bw()
  } else {
    p <- p + theme
  }

  if(do.hover){
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    return(p)
  }
}

################################# DBPlot ####################################
#' Plots continuous data on the y-axis, grouped in a settable way in the x-axis
#'
#' @param var                    Target Variable = values, OR a gene or metadata in "quotes". REQUIRED.
#' @param object                 the Seurat or RNAseq Object = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by               "metadata" to use for separating values. REQUIRED.
#' @param color.by               "metadata" to use for coloring. Affects boxplot and vlnplot fills. REQUIRED when using either.
#' @param shape.by               "metadata" to use for setting the shape of jitter.  Default = just dots. Ignored if not a "quoted" metadata or "ident"
#' @param cells.use              Cells to include: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. USE in object@cell.names[USE])
#' @param plots                  types of plots to include: possibilities = "jitter", "boxplot", "vlnplot". NOTE: The order matters, so use c("back","middle","front") when inputing multiple to put them in the order you want.
#' @param data.type              For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param do.hover               TRUE/FALSE. Default = F.  If set to true (and if there is a jitter plotted - the data it will work with) : object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  'data.hover' argument is used to determine what data to use.
#' @param data.hover             list of variable names, c("meta1","gene1","meta2","gene2"). determines what data to show on hover when do.hover is set to TRUE.
#' @param color.panel            the set of colors to draw from
#' @param colors                 indexes / or order, of colors from color.panel to actual use
#' @param main                   plot title
#' @param sub                    plot subtitle
#' @param theme                  Allows setting of a theme. Default = theme_classic when nothing is provided.
#' @param ylab                   y axis label, default is "var" or "var expression" if var is a gene.
#' @param y.breaks               a list of breaks that should be used as major gridlines. c(break1,break2,break3,etc.) NOTE: The low and highs of this variable will override `min` and `max`.
#' @param min                  Use to set a custom minimum y-value to show.  Default = set based on the limits of the data in var.
#' @param max                  Use to set a custom maximum y-value to show.  Default = set based on the limits of the data in var.
#' @param xlab                   x axis label, default is blank (NULL)
#' @param labels                 c("label1","label2","label3"). = Names of the samples/groups if you wish to change them.  Default is the values in the group.by data. NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param rotate.labels          TRUE/FALSE. whether the labels should be rotated.  Default = FALSE = vertical labels.
#' @param hline                  y value(s) where a dashed horizontal line should go
#' @param hline.linetype         Type of line.  Any ggplot linetype should work.  Defaults to "dashed"
#' @param hline.color            color(s) of the horizontal line(s)
#' @param jitter.size            the size of the jitter shapes.
#' @param jitter.width           the width/spread of the jitter in the x direction
#' @param jitter.color           the color of the jitter shapes
#' @param jitter.shapes          the shapes to use.  Default is the first in the list, which corresponds to dots.
#' @param jitter.shape.legend.size     Changes the size of the shape key in the legend.  Use a number.  OR, set to "none" to remove from the legend completely
#' @param boxplot.width          the width/spread of the boxplot in the x direction
#' @param boxplot.color          the color of the lines of the boxplot
#' @param boxplot.show.outliers  whether outliers should by including in the boxplot. Default is FALSE when there is a jitter plotted, TRUE if no jitter.
#' @param boxplot.fill          whether the boxplot should be filled in or not.
#' @param reorder.x              sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order for the samples.  Default = alphabetical then numerical.
#' @param title.legend           whether to leave the title for the plot's legend
#' @param auto.title              TRUE/FALSE = whether a default 'main' should be generated.
#' @return Makes a plot where continuous data, grouped by sample, age, cluster, etc., on the x-axis is shown on the y-axis by a violin plot, boxplot, and/or dots (or other shapes)
#' @examples
#' pbmc <- Seurat::pbmc_small
#' DBPlot("CD14", object = "pbmc", group.by = "res.1", color.by = "res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBPlot("CD14", group.by = "res.1", color.by = "res.1")

DBPlot <- function(var, object = DEFAULT, group.by, color.by,
                   shape.by = "", cells.use = NULL, plots = c("jitter","vlnplot"),
                   data.type = "normalized", do.hover = FALSE, data.hover = var,
                   color.panel = MYcolors, colors = c(1:length(color.panel)),
                   theme = theme_classic(), main = NULL, sub = NULL,
                   ylab = "make", y.breaks = NULL, min = NULL, max = NULL,
                   xlab = NULL, labels = NULL, rotate.labels = TRUE,
                   hline=NULL, hline.linetype = "dashed", hline.color = "black",
                   jitter.size=1, jitter.width=0.2, jitter.color = "black", jitter.shapes=c(16,15,17,23,25,8),
                   jitter.shape.legend.size = 3,
                   boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = NA, boxplot.fill =T,
                   reorder.x = 1:length(meta.levels(group.by, object)),
                   title.legend = F, auto.title=T){

  #Turn the object into a "name" if a full object was given
  object <- S4_2string(object)

  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  ###Determine what the y-axis should be.
  #non-direct input options are: the name of a metadata or a gene (in "quotes").
  # Both these options also get a default axis title.
  Y <- var_OR_get_meta_or_gene(var, object, data.type)

  #Unless ylab has been changed from default "make", set default y-labels if var is "metadata" or "gene".
  # If set to "var" use the 'var'.  If set to 'NULL', ylab will be removed later.
  if (!(is.null(ylab))){
    if (ylab=="make" & length(var)==1) {
      if(is.meta(var, object)){ylab <- var}
      if(is.gene(var,object)){ylab <- paste0(var," expression")}
    }
    if(ylab=="var") {ylab <- var}
    if(ylab=="make") {ylab <- NULL} #If it has not been "made" set to NULL
  }
  #Update a holder 'Shape' variable if 'shape.by' is a meta
  if (is.meta(shape.by, object)){Shape <- meta(shape.by, object)} else {Shape = NA}

  #Groundwork for plotly hover data:
  #Overall: if do.hover=T and data.hover has a list of genes / metas,
  # then for all cells, make a string "var1: var1-value\nvar2: var2-value..."
  hover.string <- NA
  if (do.hover) {
    hover.string <- make_hover_strings(data.hover, object, data.type)
  }

  ###Make dataframe for storing the plotting data:
  #The data =Y, how to group the data = sample, how to color the groupings = color, and what the shape should be if there is a "jitter" made.
  full_dat <- data.frame(Y = Y,
                         sample = meta(group.by, object),
                         color = meta(color.by, object),
                         shape = Shape,
                         hover.string = hover.string)
  #Subset the data.frame to only the cells in cell.use.
  Target_dat <- full_dat[all.cells %in% cells.use,]

  #Reorder x groupings (steps 1 and 2)
  #1-Rename the grouping (=Target_dat$sample) labels in order to set their order.
  #2-Store originals in orig.names.
  #3-Names will be set labels or back to orig.names further down in the code.
  if (typeof(reorder.x)=="integer" | typeof(reorder.x)=="double"){
    #This step is necessary becuase ggplot orders by "character" which would put 10 after 1 instead of 9.
    reorder.x <- as.character(unlist(sapply(reorder.x, function(X) ifelse(X<10,paste0("0",X),X))))
  }
  orig.names <- levels(as.factor(as.character(Target_dat$sample)))
  Target_dat$sample <- as.factor(as.character(Target_dat$sample))
  levels(Target_dat$sample) <- reorder.x
  Target_dat$sample <- as.factor(as.character(Target_dat$sample))

  #####Start making the plot
  p <- ggplot(Target_dat, aes(x=sample, y=Y, fill=color)) + theme
  #Set the legend to not have a title, unless requested.
  if (!title.legend) {p <- p + theme(legend.title=element_blank())}
  #Add the y label
  p<- p + ylab(ylab)
  #Set the y-axis limits if a min or max is given.
  if ((!(is.null(min))) & (!(is.null(max)))){
    p <- p + ylim(min,max)
  } else {
    if (!(is.null(min))){
      p <- p + ylim(min, max(Target_dat$Y))
    }
    if (!(is.null(max))){
      p <- p + ylim(min(Target_dat$Y), max)
    }
  }

  ###Add data based on what is requested in plots, *ordered by their order*
  for (i in 1:length(plots)){
    #If next request is "boxplot", make a boxplot.
    if (plots[i] == "boxplot") {
      if (is.na(boxplot.show.outliers)){
        boxplot.show.outliers <- ifelse("jitter" %in% plots, FALSE, TRUE)
      }
      if (boxplot.show.outliers) {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0))
      } else {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0),
                              outlier.shape = NA)
      }
    }
    #If next request is "jitter", make a jitter.  a.k.a. add dots with y=Y, and randomized spread in x direction.
    if (plots[i] == "jitter") {
      #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
      if (is.meta(shape.by, object)){
        #Make jitter with shapes
        p <- p + geom_jitter(size=jitter.size,
                             width=jitter.width,
                             height = 0,
                             #Add the hover info
                             if(do.hover){aes(shape = shape, text = hover.string)}
                             else{aes(shape = shape)},
                             color = jitter.color)
        #Actually set the shapes to jitter.shapes.
        p <- p + scale_shape_manual(values = jitter.shapes[1:length(levels(as.factor(Target_dat$shape)))],
                                    labels = levels(as.factor(as.character(Target_dat$shape))))
        #Also change the size of the shape key in the legend, unless jitter.shape.legend.size was set to NA or "none".
        if (!is.na(jitter.shape.legend.size) & jitter.shape.legend.size!="none"){
          p <- p + guides(shape = guide_legend(override.aes = list(size=jitter.shape.legend.size)))
        }
        if (jitter.shape.legend.size=="none"){
          p <- p + guides(shape = "none")
        }
      } else {
        p <- p + geom_jitter(size=jitter.size,
                             width=jitter.width,
                             height = 0,
                             if(do.hover){aes(text = hover.string)},
                             shape = jitter.shapes[1],
                             color = jitter.color)
      }
    }
    if (plots[i] == "vlnplot") {
      p <- p + geom_violin()
    }
  }

  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype= hline.linetype, color = hline.color)}

  #Set default title if no 'main' was given and if 'autotitle' is set to TRUE
  if (is.null(main) & auto.title==T & length(var)==1){
    if(is.meta(var, object)){main <- var}
    if(is.gene(var, object)){main <- var}
  }
  #Change the x-axis labels if a labels was given.
  if (!(is.null(labels))){
    p <- p + scale_x_discrete(labels=labels)
  } else {
    # If it was not given, put the orig.names back
    p <- p + scale_x_discrete(labels=orig.names[order(reorder.x)])
  }
  #Change y-axis limits/breaks if requested
  if (!is.null(y.breaks)) {
    p <- p + scale_y_continuous(breaks = y.breaks) + coord_cartesian(ylim=c(min(y.breaks),max(y.breaks)))
  }
  #Rotate Labels if rotate.labels = TRUE
  if (rotate.labels) {p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1, size=12))}
  #Set colors to color.panel[colors], set the x-axis formatting, and add titles.
  p <- p + scale_fill_manual(values=color.panel[colors]) +
    ggtitle(main, sub) + xlab(xlab)

  #DONE. Return the plot
  if(do.hover & ("jitter" %in% plots)){
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    return(p)
  }
}

########## DBBarPlot: Builds a stacked bar plot to show the composition of samples / ages / 'group.by' ##########
#' Outputs a stacked bar plot to show the percent composition of samples or other cell groupings
#'
#' @param var                    Target Variable = values, OR a metadata in "quotes". REQUIRED. Length must be same  as the number of cells in the Seurat object or Samples in the RNAseq object
#' @param object                 the Seurat or RNAseq object to draw from = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by               "metadata" to use for separating values. REQUIRED (Default is to use a metadata named "Sample").
#' @param cells.use              Cells to include: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. *USE* in object@cell.names[*USE*])
#' @param color.panel            the set of colors to draw from
#' @param colors                 indexes / or order, of colors from color.panel to actual use
#' @param do.hover               TRUE/FALSE. Default = F.  If set to true, object will be converted to a ggplotly object so that data about individual bars will be displayed when you hover your cursor over them.  Data displayed will be the "counts" and "percentage of total".'data.hover' argument is not used with this plotting function.
#' @param theme                  Allows setting of a theme. Default = theme_classic() when nothing is provided.
#' @param xlab                   "character". The text title for the x axis.  NULL/blank by default.
#' @param ylab                   "character". The text title for the y axis.  Auto-generated by default.  Provide ylab = NULL to remove
#' @param x.labels               Replacement x-axis labels to use instead of the identities of gorup.by
#' @param rotate.labels          TRUE/FALSE. whether the labels should be rotated.  Default = FALSE = vertical labels.
#' @param y.breaks               The numerical labels for the y axis.  Note: The percentages that build this figure will always add up to 1, so the plot will always go from 0 to 1.
#' @param main                   "character". Plot main title.
#' @param sub                    "character". Plot subtitle.
#' @param rename.groups          new names for the identities of var.  Change to NULL to remove labeling altogether.
#' @param legend.title           Title for the legend.  Default = blank / NULL
#' @param reorder.x              sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order for the samples.  Default = alphabetical then numerical.
#' @return Makes a plot where discrete data is shown on the y-axis as "percent of total" in a stacked barplot, grouped by sample, age, cluster, etc., on the x-axis.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' DBBarPlot("ident", object = "pbmc", group.by = "orig.ident")
#' #For real data, you will have more distinct samples than this truncated dataset.
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBBarPlot("ident", group.by = "orig.ident")

DBBarPlot <- function(var="ident", object = DEFAULT, group.by = "Sample",
                      cells.use = NULL,
                      color.panel = MYcolors, colors = c(1:length(color.panel)),
                      do.hover = F, theme = theme_classic(),
                      xlab = NULL, ylab = "make", x.labels = NA, rotate.labels = TRUE,
                      y.breaks = c(0,0.5,1),
                      main = "make", sub = NULL, rename.groups = NA,
                      legend.title = NULL,
                      reorder.x = 1:length(meta.levels(group.by, object))
                      ){

  #Turn the object into a "name" if a full object was given
  object <- S4_2string(object)

  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  ####Retrieve metas: var, group.by
  #var to y.var
  #If name of meta in "quotes", obtain the meta
  if(length(var)==1 & typeof(var)=="character") {
    if (is.meta(var, object)){
      y.var <- as.factor(meta(var, object))
    }
  } else {y.var <- var}
  #group.by to x.var
  #If name of meta in "quotes", obtain the meta
  if(length(group.by)==1 & typeof(group.by)=="character") {
    if (is.meta(group.by, object)){
      x.var <- as.factor(meta(group.by, object))
    }
  }
  #Subset the x.var and y.var to only the cells in cell.use.
  x.var <- as.factor(as.character(x.var[all.cells %in% cells.use]))
  y.var <- y.var[all.cells %in% cells.use]
  #Reorder x groupings (steps 1 and 2)
  #1-Rename the x.var labels in order to set their order.
  #2-Store originals in orig.names.
  #3-Names will be set back to orig.names or x.labels further down on in the code.
  if (typeof(reorder.x)=="integer"){
    reorder.x <- as.character(unlist(sapply(reorder.x, function(X) ifelse(X<10,paste0("0",X),X))))
  }
  orig.names <- levels(x.var)
  levels(x.var) <- reorder.x

  #Groundwork for plotly hover data:
  #Overall: if do.hover=T and data.hover has a list of genes / metas,
  # then for all cells, make a string "var1: var1-value\nvar2: var2-value..."
  hover.string <- NA
  if (do.hover) {
    features.info <- data.frame(name = rep(meta.levels(var), length(levels(x.var))),
                                y.counts = c(sapply(levels(as.factor(x.var)), function(X)
                                  unlist(sapply(levels(as.factor(y.var)), function(Y)
                                    #Number of Xs that are Ys
                                    sum(y.var==Y & x.var == X))))),
                                y.percents = c(sapply(levels(as.factor(x.var)), function(X)
                                  unlist(sapply(levels(as.factor(y.var)), function(Y)
                                    #Number of Xs that are Ys, divided by the total number of Xs.
                                    sum(y.var==Y & x.var == X)/sum(x.var == X))))))
    names(features.info)<-c("Identity","Count","Percent of total")
    hover.string <- sapply(1:nrow(features.info), function(row){
      paste(as.character(sapply(1:ncol(features.info), function(col){
        paste0(names(features.info)[col],": ",features.info[row,col])})),collapse = "\n")
    })
  }

  #Build data (Make a dataframe while calculating the percent makeup of x.var groups by y.var identities.)
  #Generate the x.grouping data (needs to be the identities of x.var each individually repeated
  # the number of times that there are distinct levels in the var / y.var.)
  dat <- data.frame(grouping = c(sapply(levels(as.factor(x.var)), function(X) rep(X, length(levels(as.factor(y.var)))))),
                    #Label the y.ident that this perentage is for. (Used for coloring)
                    y.ident = rep(levels(as.factor(y.var)), length(levels(as.factor(x.var)))),
                    #Generate the percents
                    y.percents = c(sapply(levels(as.factor(x.var)), function(X)
                      unlist(sapply(levels(as.factor(y.var)), function(Y)
                        #Number of Xs that are Ys, divided by the total number of Xs.
                        sum(y.var==Y & x.var == X)/sum(x.var == X)
                      ))
                    )),
                    hover.string = hover.string
  )

  #Build Plot
  p <- ggplot(data=dat, aes(x = grouping)) + theme +
    #Add the bars.
    if(do.hover){
      geom_col(aes(y=y.percents, fill = y.ident, text = hover.string))
    } else {
      geom_col(aes(y=y.percents, fill = y.ident))
    }
  #Populate ylab if left as "make".
  if(ylab == "make"){ ylab <- paste0("Percent of ",
                                     ifelse(grepl("Seurat",classof(object)),
                                            "cells",
                                            "samples"))
  }
  #Add the y-axis labels and name to the plot
  p <- p + scale_y_continuous(breaks= y.breaks,
                              limits = c(0,1),
                              name = ylab)
  ###Set the colors & and rename the groupings
  #If labels.rename was changed from NA, rename the labels to rename.groups
  if (!(is.na(rename.groups[1]))){
    p <- p+ scale_fill_manual(name = legend.title,
                              values = color.panel[colors],
                              labels = rename.groups)
  } else {
    #If not, just set the colors and name for the legend key
    p <- p+ scale_fill_manual(name = legend.title,
                              values = color.panel[colors])
  }

  #Set var to default if 'main' left as "make".
  if (!(is.null(main))){
    if (main=="make" & length(var)==1){
      main <- var
    }
    if (main=="make"){
      main <- NULL
    }
  }

  #Add the x axis title, plot title and subtitle,
  p <- p + xlab(xlab) +
    ggtitle(main, subtitle = sub) +
    theme(axis.text.x= element_text(size=10)) +
    if (rotate.labels) {theme(axis.text.x= element_text(angle=45, hjust = 1.3, vjust = 1.2, size=12))} +
    theme(legend.title=element_text(size=12)) +
    theme(legend.text=element_text(size=14))

  #Rename the x-axis labels if an x.labels was given.
  if (!(is.na(x.labels[1]))){
    p <- p + scale_x_discrete(labels=x.labels)
  } else {
    # If it was not given, put the orig.names back
    p <- p + scale_x_discrete(labels=orig.names[order(reorder.x)])
  }

  #DONE. Return the plot
  if(do.hover){
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    return(p)
  }
}

########## DBHeatmap: Builds a heatmap with given genes using a pretty fast heatmap algorithm.  Unfortunately, it's a bit buggy. ##########
#' Outputs a heatmap of the object
#'
#' @param genes c("gene1","gene2","gene3",...) = list of genes to put in the heatmap. REQUIRED.
#' @param object the Seurat or RNAseq object to draw from = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param cells.use Cells to include: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. *USE* in object@cell.names[*USE*])
#' @param data.type For obtaining the expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param cell.names.meta quoted "name" of a meta.data slot to use for naming the columns instead of the raw cell/sample names.  Default = just use the raw cell/sample names in the data slot.
#' @param heatmap.colors the colors to use for the expression matrix. Default is a ramp from navy to white to red with 50 slices.
#' @param cells.annotation FALSE/ the name of a metadata slot containing how the cells/samples should be annotated. Default = FALSE.
#' @param cells.annotation.colors The colors to use for the annotations bar
#' @param data.out = If set to TRUE, the output will be a list of objects and and the script that would have been used for inputing into pheatmap.  NOTE: currently, a row#### arguement is added that should not be there.  I will remove this in a later release, but for now, know that you should remove this before trying to run the script.
#' @param ... other arguments that will be passed to pheatmap. For scRNAseq data, I recommend running with cluster_cols=FALSE first because clustering of thousands of cells can tak a long time!  Common ones: main for adding a title, cluster_cols or cluster_rows = FALSE to turn clustering off for cells or genes, treeheight_row or treeheight_col = 0 to remove or a positive # for setting how large the trees on the side/top should be drawn
#' @description This function is a wrapper for the pheatmap function of the pheatmap package.  Given a set of genes, it will spit out a basic heatmap from your RNAseq data, or just certain cells/samples.
#' @return Makes a heatmap where
#' @examples
#' pbmc <- Seurat::pbmc_small
#' DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'             "FCGR3A","LYZ","PPBP","CD8A"),
#'           object = "pbmc",
#'           cells.annotation = "ident")
#' #For real data, you will have more cells than this truncated dataset,
#' # so I recommend turning off cell clustering when you are trying out tweaks by
#' # adding cluster_cols=FALSE
#' DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'             "FCGR3A","LYZ","PPBP","CD8A"),
#'           object = "pbmc",
#'           cells.annotation = "ident",
#'           cluster_cols=FALSE)
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
#'             "FCGR3A","LYZ","PPBP","CD8A"),
#'           cells.annotation = "ident")
#'
DBHeatmap <- function(genes=NULL, object = DEFAULT, cells.use = NULL,
                      cell.names.meta = NULL, data.type = "normalized",
                      heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50),
                      cells.annotation = FALSE, cells.annotation.colors = rep(list(MYcolors),length(cells.annotation)),
                      data.out=FALSE, ...){

  #If no genes given, "error out"
  if(is.null(genes)){
    return('This function is not yet set up to select which genes to use.\nPlease provide a list of genes as genes=c("gene1","gene2","gene3")')
  }

  #Turn the object into a "name" if a full object was given
  object <- S4_2string(object)

  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  #Make the data matrix
  data <- as.matrix(which_data(data.type,object)[genes,cells.use])

  #Make the cells.annotations data for color annotation bar
  Col_annot <- NA
  Col_annot_colors <- NA
  if(!(cells.annotation[1] == FALSE)){
    Col_annot <- data.frame(as.factor(as.character(meta(cells.annotation[1], object)[all.cells %in% cells.use])),
                            row.names = cells.use)
    names(Col_annot) <- cells.annotation[1]
    new_cols <- unlist(cells.annotation.colors[[1]])[1:length(levels(Col_annot[[1]]))]
    names(new_cols) <- levels(Col_annot[[1]])
    Col_annot_colors <- list(new_cols)
    names(Col_annot_colors) <- cells.annotation[1]
    if (length(cells.annotation)>1){
      print(paste0("Currently, only one annotation per heatmap is supported. To create a plot with two or more, you can create separate heatmaps, and stitch them together manually in other software. Generating heatmap with '", cells.annotation[1], "'..."))
      # for (i in 2:length(cells.annotation[-1])){
      #   eval(expr = parse(text = paste0("Col_annot$'",cells.annotation[i],"' <- as.factor(as.character(meta(cells.annotation[i],object)[all.cells %in% cells.use]))")))
      #   new_cols <- unlist(cells.annotation.colors[[i]])[1:length(levels(Col_annot[[i]]))]
      #   names(new_cols) <- levels(Col_annot[[i]])
      #   eval(expr = parse(text = paste0("Col_annot_colors$'",cells.annotation[i],"' <- list(new_cols)")))
      # }
    }
    rownames(Col_annot) <- colnames(data)
  }

  #Establish cell/sample names
  Col_names <- NULL
  if(!(is.null(cell.names.meta))){
    Col_names <- as.character(meta(cell.names.meta, object)[all.cells %in% cells.use])
  }

  if(data.out){
    OUT <- list(data,heatmap.colors,Col_annot,Col_annot_colors, Col_names,
                paste0("pheatmap(mat = data, color = heatmap.colors, annotation_col = Col_annot, ",
                       "annotation_colors = Col_annot_colors, scale = 'row', ",
                       "labels_col = Col_names"))
    names(OUT)<- c("data","heatmap.colors","Col_annot","Col_annot_colors","Col_names","pheatmap.script")
    return(OUT)
  }

  #Make the Heatmap
  HM <- pheatmap::pheatmap(mat = data,
                           color = heatmap.colors,
                           annotation_col = Col_annot,
                           annotation_colors = Col_annot_colors,
                           scale = "row",
                           labels_col = Col_names,
                           ...)

  #DONE: Return the heatmap
  return(HM)
}

#### DBPlot_multi_var_summary : Generates a DBPlot where datapoints are genes/metadata instead of cells/samples.
#' Generates a DBPlot where datapoints are gene(s)/metadata instead of cells/samples.
#'
#' @param vars               c("gene1","gene2","gene3",...). REQUIRED. A list of genes from which to generate the plot
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by           "metadata" to use for separating values. REQUIRED.
#' @param color.by               "metadata" to use for coloring. Affects boxplot and vlnplot fills. REQUIRED when using either.
#' @param cells.use              Cells to include in generating the summary data: either in the form of a character list of names, or a logical that is the same length as the number of cells in the object (a.k.a. USE in object@cell.names[USE])
#' @param plots                  types of plots to include: possibilities = "jitter", "boxplot", "vlnplot". NOTE: The order matters, so use c("back","middle","front") when inputing multiple to put them in the order you want.
#' @param data.type              In what format should the data for each var be summarized? DEFAULT = "relative".  Most other options are not as great for comparing accross genes that may have vastly different expression patterns. All options: "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)?
#' @param data.summary           "mean" or "median". = the summary statistic to use for summarizing expression/score accross the gorups. Default is to use the mean.
#' @param do.hover               TRUE/FALSE. Default = F.  If set to true: object will be converted to a ggplotly object so that data about individual points will be displayed when you hover your cursor over them.  Data shown will be the gene name. ('data.hover' argument is NOT used with this function.)
#' @param color.panel            the set of colors to draw from
#' @param colors                 indexes / or order, of colors from color.panel to actual use
#' @param theme                  Allows setting of a theme. Default = theme_classic when nothing is provided.
#' @param legend.show            TRUE/FALSE, whether the legend should be included/removed. Default = FALSE
#' @param legend.title           TRUE/FALSE/"custom-title", whether to have title for the plot's legend OR the custom title you would like the legend to have. Default = FALSE
#' @param main                   plot title (Default is none)
#' @param sub                    plot subtitle (Default is none)
#' @param ylab                   y axis label (Default is none)
#' @param y.breaks               a list of breaks that should be used as major gridlines. c(break1,break2,break3,etc.)
#' @param min                  Use to set a custom minimum y-value to show.  Default = set based on the limits of the data in var.
#' @param max                  Use to set a custom maximum y-value to show.  Default = set based on the limits of the data in var.
#' @param xlab                   x axis label, default is blank (NULL)
#' @param labels                 c("label1","label2","label3"). = Names of the samples/groups if you wish to change them.  Default is the values in the group.by data. NOTE: you need to give at least as many labels as there are discrete values in the group.by data.
#' @param rotate.labels          TRUE/FALSE. whether the labels should be rotated.  Default = FALSE = vertical labels.
#' @param reorder.x              Sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order for the samples.  Default = alphabetical then numerical. Method: look at an un-reordered plot, then plug in which position from left-to-right that you want the left-most group to move to, then where the originally second-from-left should go, and so on.
#' @param jitter.size            the size of the jitter shapes.
#' @param jitter.width           the width/spread of the jitter in the x direction
#' @param jitter.color           the color of the jitter shapes
#' @param jitter.shape           number corresponding to a ggplot shape.  Default = 16 = a small dot.
#' @param boxplot.width          the width/spread of the boxplot in the x direction
#' @param boxplot.color          the color of the lines of the boxplot
#' @param boxplot.show.outliers  whether outliers should by including in the boxplot. Default is FALSE when there is a jitter plotted, TRUE if no jitter.
#' @param boxplot.fill          whether the boxplot should be filled in or not.
#' @param hline                  y value(s) where a dashed horizontal line should go
#' @param hline.linetype         Type of line.  Any ggplot linetype should work.  Defaults to "dashed"
#' @param hline.color            color(s) of the horizontal line(s)
#' @return This function will output a DBPlot where each data point represents the average (or median) of expression of an individual gene, or of an individual metadata score, across a group of samples.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' genes <- get.genes("pbmc")[1:30]
#' DBPlot_multi_var_summary(genes, object = "pbmc",
#'                          group.by = "res.1", color.by = "res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object
#'       # input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBPlot_multi_var_summary(genes,
#'                          group.by = "res.1", color.by = "res.1")
#'
#' # To change it to have the violin plot in the back, a jitter on
#' #  top of that, and a white boxplot with no fill in front:
#' DBPlot_multi_var_summary(genes, object = "pbmc",
#'                          group.by = "res.1", color.by = "res.1",
#'                          plots = c("vlnplot","jitter","boxplot"),
#'                          boxplot.color = "white", boxplot.fill = F)

DBPlot_multi_var_summary <- function(vars, object = DEFAULT, group.by="Sample", color.by=NULL, cells.use = NULL,
                                     plots = c("vlnplot","jitter"), data.type = "relative", data.summary = "mean",
                                     do.hover = FALSE,
                                     color.panel = MYcolors, colors = c(1:length(color.panel)),
                                     theme = theme_classic(), legend.show = FALSE, legend.title = FALSE,
                                     main = NULL, sub = NULL,
                                     ylab = NULL, y.breaks = NULL, min = NULL, max = NULL,
                                     xlab = NULL, labels = NULL, rotate.labels = TRUE,
                                     reorder.x = 1:length(meta.levels(group.by, object)),
                                     jitter.size=1, jitter.width=0.2, jitter.color = "black",
                                     jitter.shape = 16,
                                     boxplot.width = 0.2, boxplot.color = "black",
                                     boxplot.show.outliers = NA, boxplot.fill =T,
                                     hline=NULL, hline.linetype = "dashed", hline.color = "black"){

  #Turn the object into a "name" if a full object was given
  object <- S4_2string(object)

  #Populate cells.use with a list of names if it was given anything else.
  cells.use <- which_cells(cells.use, object)
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)

  #### Ensure that vars is a list of numerical values.

  #### Run prep for setting color
  if(is.null(color.by)){
    color.by <- group.by
  }
  #### Ensure that there are no ambiguities between group.by and color.by
  grp.color.check.matrix <- as.matrix(table(meta(group.by,object),meta(color.by,object)))
  grp.color.check.matrix <- data.frame(grp.color.check.matrix>0)
  if(sum(rowSums(grp.color.check.matrix)!=1)){return(print("Unable to interpret color.by input.  group.by and color.by metadatas must match so that only one color will be assigned to each group."))}
  color <- sapply(seq_len(dim(grp.color.check.matrix)[1]), function(X)
    names(grp.color.check.matrix)[grp.color.check.matrix[X,]==TRUE])
  names(color) <- row.names(grp.color.check.matrix)


  #### Create data table
  # Summarizes data and creates vars x groupings table
  groupings <- as.factor(meta(group.by, object)[cells.use %in% all.cells])
  if(data.summary=="mean"){
    summarys <- data.frame(sapply(levels(groupings), function(this.group)
      sapply(vars, function(X) mean(var_OR_get_meta_or_gene(X,object,data.type)[groupings==this.group])
      )))
  } else {
    if(data.summary=="median"){
      summarys <- data.frame(sapply(levels(groupings), function(this.group)
        sapply(vars, function(X) median(var_OR_get_meta_or_gene(X,object,data.type)[groupings==this.group])
        )))
    } else {
      return(print("mean and median are the only summary statistics currently supported.", quote = F))
    }
  }

  # Turn that data.table into a linear list that is part of a ggplot friendly data.frame
  dat <- data.frame(Y = unlist(summarys),
                    vars = rep(rownames(summarys),ncol(summarys)),
                    X = c(sapply(names(summarys), function(X) rep(X, nrow(summarys)))))
  dat$color <- color[dat$X]

  #Reorder x groupings (steps 1 and 2, #3 is completed later in the code)
  if (is.numeric(reorder.x)){
    #This step is necessary becuase ggplot orders by "character" which would put 10 after 1 instead of 9.
    reorder.x <- as.character(unlist(sapply(reorder.x, function(X) ifelse(X<10,paste0("0",X),X))))
  }
  #1-Store originals in orig.names.
  orig.names <- levels(as.factor(as.character(dat$X)))
  #2-Rename the groupings (=dat$X) in order to set their order.
  dat$X <- as.factor(as.character(dat$X))
  levels(dat$X) <- reorder.x
  dat$X <- as.factor(as.character(dat$X))
  #3-Names will be set back to orig.names (unless new labels were provided) further down in the code.

  ####If do.hover = T, make the string for implementing it
  #Make data.frame genes x display data
  features.info <- data.frame(gene = dat$vars)
  hover.string <- sapply(seq_len(nrow(features.info)), function(row){
    paste(as.character(sapply(seq_len(ncol(features.info)), function(col){
      paste0(names(features.info)[col],": ",features.info[row,col])})),collapse = "\n")
  })

  #####Start making the plot
  p <- ggplot(dat, aes(x=X, y=Y, fill=color)) + theme
  #Set the y-axis limits if a min or max is given.
  if ((!(is.null(min))) & (!(is.null(max)))){
    p <- p + ylim(min,max)
  } else {
    if (!(is.null(min))){
      p <- p + ylim(min, max(Target_dat$Y))
    }
    if (!(is.null(max))){
      p <- p + ylim(min(Target_dat$Y), max)
    }
  }
  #Add data based on what is requested in plots, *ordered by their order*
  for (i in 1:length(plots)){
    #If next request is "boxplot", make a boxplot.
    if (plots[i] == "boxplot") {
      if (is.na(boxplot.show.outliers)){
        boxplot.show.outliers <- ifelse("jitter" %in% plots, FALSE, TRUE)
      }
      if (boxplot.show.outliers) {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0))
      } else {
        p <- p + geom_boxplot(width=boxplot.width, color = boxplot.color,
                              alpha = ifelse(boxplot.fill, 1, 0),
                              outlier.shape = NA)
      }
    }
    #If next request is "jitter", make a jitter.  a.k.a. add dots with y=Y, and randomized spread in x direction.
    if (plots[i] == "jitter") {
      p <- p + geom_jitter(size=jitter.size,
                           width=jitter.width,
                           height = 0,
                           if(do.hover){aes(text = hover.string)},
                           shape = jitter.shape,
                           color = jitter.color)
    }
    if (plots[i] == "vlnplot") {
      p <- p + geom_violin()
    }
  }
  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype= hline.linetype, color = hline.color)}
  #Change y-axis limits/breaks if requested
  if (!is.null(y.breaks)) {
    p <- p + scale_y_continuous(breaks = y.breaks)
  }
  #Rotate Labels if rotate.labels = TRUE
  if (rotate.labels) {p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1, size=12))}
  #Set colors to color.panel[colors], and add/remove legend title.
  if (is.logical(legend.title)){
    if(!legend.title) {
      #legend.title = FALSE
      p <- p + scale_fill_manual(values=color.panel[colors]) + theme(legend.title=element_blank())
    } else {
      #legend.title = TRUE
      p <- p + scale_fill_manual(name = color.by,
                                 values = color.panel[colors])
    }
  } else {
    #legend.title = "custom"
    p <- p + scale_fill_manual(name = legend.title,
                               values = color.panel[colors])
  }
  #Change the x-axis labels if a labels was given.
  if (!(is.null(labels))){
    p <- p + scale_x_discrete(labels=labels)
  } else {
    #If new labels were not given, put the original names back.
    p <- p + scale_x_discrete(labels=orig.names[order(reorder.x)])
  }
  #Set titles and x-axis label
  p <- p + ggtitle(main, sub) + xlab(xlab) + ylab(ylab)
  #Remove legend unless wanted
  p <- p + theme(legend.position = ifelse(legend.show, "right", "none"))

  #DONE. Return the plot
  if(do.hover & ("jitter" %in% plots)){
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    return(p)
  }
}

#### multiDBPlot : a function for quickly making multiple DBPlots arranged in a grid.
#' Generates multiple DBPlots arranged into a grid.
#'
#' @param vars               c("var1","var2","var3",...). REQUIRED. A list of vars from which to generate the separate plots
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by           "metadata" to use for separating values. REQUIRED.
#' @param color.by           "metadata" to use for coloring. Affects boxplot and vlnplot fills. REQUIRED when using either.
#' @param show.legend        TRUE/FALSE. Whether or not you would like a legend to be plotted.  Default = FALSE
#' @param ncol               #. How many plots should be arranged per row.  Default = 3.
#' @param nrow               #/NULL. How many rows to arrange the plots into.  Default = NULL(/blank) --> becomes however many rows are needed to show all the data.
#' @param add.title          TRUE/FALSE. Whether a title should be added.
#' @param ylab               TRUE/FALSE/"var". Whether a y-axis title should be added, TRUE/FALSE or "var".  If "var", then the var names will be used.  If TRUE, y-label for any gene vars will be "'var' expression"
#' @param OUT.List           TRUE/FALSE. (Default = FALSE) Whether the output should be a list of objects instead of the full plot.  Outputting as list allows manual input into gridArrange for moving plots around / adjusting sizes.  In the list, all plots will be named by the variable being shown.
#' @param ...                other paramters that can be given to DBPlot function used in exactly the same way.
#' @return Given multiple 'var' parameters, this function will output a DBPlot for each one, arranged into a grid.  All parameters that can be adjusted in DBPlot can be adjusted here.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' genes <- c("CD8A","CD3E","FCER1A","CD14")
#' multiDBPlot(genes, object = "pbmc",
#'             group.by = "res.1", color.by = "res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object
#'       # input can be skipped completely.
#' DEFAULT <- "pbmc"
#' multiDBPlot(genes,
#'             group.by = "res.1", color.by = "res.1")
#' #To make it output a grid that is 2x2, to add y-axis labels
#' # instead of titles, and to show legends...
#' multiDBPlot(genes,
#'             group.by = "res.1", color.by = "res.1",
#'             nrow = 2, ncol = 2,              #Make it 2x2
#'             add.title = FALSE, ylab = TRUE,  #Add y axis labels instead of titles
#'             show.legend = TRUE)              #Show legends
#' # To eliminate the "expression", change ylab = TRUE to ylab = "var"
#' multiDBPlot(genes,
#'             group.by = "res.1", color.by = "res.1",
#'             nrow = 2, ncol = 2,              #Make it 2x2
#'             add.title = FALSE, ylab = "var", #Add y axis labels without "expression"
#'             show.legend = TRUE)              #Show legends

multiDBPlot <- function(vars, object = DEFAULT, group.by, color.by,
                        show.legend = FALSE,
                        ncol = 3,
                        nrow = NULL,
                        add.title=TRUE,
                        ylab = FALSE,
                        OUT.List = FALSE,
                        ...){

  ylab.input <- ylab

  plots <- lapply(vars, function(X) { DBPlot(X, object, group.by, color.by,
                                             ylab = if(typeof(ylab.input)=="character"){ifelse((ylab.input == "var"),X,return('Error: ylab must be TRUE/FALSE or "var"'))} else {if(ylab.input){"make"}else{NULL}},
                                             auto.title = add.title,
                                             ...) +
      theme(legend.position = ifelse(show.legend, "right", "none"))
  })

  #Output
  if (OUT.List){
    names(plots) <- vars
    return(plots)
  } else {
    return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
  }
}

#### multiDBDimPlot : a function for quickly making multiple DBDimPlots arranged in a grid.
#' Generates multiple DBDimPlots arranged into a grid.
#'
#' @param vars               c("var1","var2","var3",...). REQUIRED. A list of vars from which to generate the separate plots
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param show.legend        TRUE/FALSE. Whether or not you would like a legend to be plotted.  Default = FALSE
#' @param ncol               #. How many plots should be arranged per row.  Default = 3.
#' @param nrow               #/NULL. How many rows to arrange the plots into.  Default = NULL(/blank) --> becomes however many rows are needed to show all the data.
#' @param add.title          TRUE/FALSE. Whether a title should be added.
#' @param axes.labels        TRUE/FALSE. Whether a axis labels should be added
#' @param OUT.List           TRUE/FALSE. (Default = FALSE) Whether the output should be a list of objects instead of the full plot.  Outputting as list allows manual input into gridArrange for moving plots around / adjusting sizes.  In the list, all plots will be named by the variable being shown.
#' @param ...                other paramters that can be given to DBDimPlot function used in exactly the same way.
#' @return Given multiple 'var' parameters, this function will output a DBDimPlot for each one, arranged into a grid.  All parameters that can be adjusted in DBDimPlot can be adjusted here, but the only parameter that can be adjusted between each is 'var'.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' genes <- c("CD8A","CD3E","FCER1A","CD14","MS4A1")
#' multiDBDimPlot(c(genes, "ident"), object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' multiDBDimPlot(c(genes, "ident"))

multiDBDimPlot <- function(vars, object = DEFAULT,
                           show.legend = F,
                           ncol = 3, nrow = NULL,
                           add.title=T, axes.labels=F,
                           OUT.List,
                           ...){

  #Interpret axes.labels: If left as FALSE, set lab to NULL so they will be removed.
  # If set to TRUE, set it to "make".
  lab <- if(!axes.labels) {NULL} else {"make"}

  plots <- lapply(vars, function(X) { DBDimPlot(X, object,
                                                auto.title = add.title,
                                                xlab = lab,
                                                ylab = lab,
                                                ...) +
      theme(legend.position = ifelse(show.legend, "right", "none"))
  })
  if (OUT.List){
    names(plots) <- vars
    return(plots)
  } else {
    return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
  }
}

##################### multiDBDimPlot_vary_cells #######################
#' Generates multiple DBDimPlots, each showing different cells, arranged into a grid.
#'
#' @param var                name of a "gene", "meta.data", "ident", or a vector the length of the #cells or #samples in 'object'. REQUIRED. A var with which to color the separate plots.  Referencing discrete or continuous data are both allowed.
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param cells.use.meta     The name of the meta.data that will be used for selecting cells. REQUIRED.
#' @param cells.use.levels   The values/groupings of the cells.use.meta that you wish to show. NOTE: these will be put into the plot's main title in order to have them be identifiable.
#' @param all.cells.plot     TRUE/FALSE, whether a plot showing all of the cells should be included at the end.
#' @param show.legend        TRUE/FALSE. Whether or not you would like a legend to be plotted in every plot.  Default = FALSE
#' @param add.single.legend  TRUE/FALSE, whether to add a single legend as an additional plot.
#' @param ncol               #. How many plots should be arranged per row.  Default = 3.
#' @param nrow               #/NULL. How many rows to arrange the plots into.  Default = NULL(/blank) --> becomes however many rows are needed to show all the data.
#' @param add.title          TRUE/FALSE. Whether a title should be added.
#' @param axes.labels        TRUE/FALSE. Whether a axis labels should be added.
#' @param main               For setting a title for each individual plot.  Options: Leave as "make" and the titles will be which cells are in each plot.  "any.text" will set the same title for all plots.  NULL or "" will remove the titles.
#' @param all.cells.main     For setting the title of the All.cells plot.  Default is All Cells.  Set to NULL or "" to remove.
#' @param min                Works as in DimPlot, not required to be given, but sets the value associated with the minimum color.  All points with a lower value than this will get the same min.color.
#' @param max                Works as in DimPlot, not required to be given, but sets the value associated with the maximum color.  All points with a higher value than this will get the same max.color.
#' @param ...                other paramters that can be given to DBDimPlot function used in the same way.
#' @param data.type          As with regular DBDimPlotting, For when plotting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @param OUT.List           TRUE/FALSE. (Default = FALSE) Whether the output should be a list of objects instead of the full plot.  Outputting as list allows manual input into gridArrange for moving plots around / adjusting sizes.  In the list, the All.cells plot will be named "all" and the legend will be named "legend". Others will be named by which cells are in them.
#' @return A function for quickly making multiple DBDimPlots arranged in a grid, where instead of varying the 'var' displayed, what varies is the cells that are shown. Most parameters that can be adjusted in DBDimPlot can be adjusted here, but the only parameter that can be adjusted between each plot is which cells get displayed. NOTE: This function is incompatible with changing the 'colors' input. If you need to change the order of when certain colrs are chosen, change the order in color.panel. Also note: if 'var' given refers to continuous data, then the min and max will be calculated at the beginning in order to make the scale consistent accross all plots.  You can still set your own range, and this would also create a consistent scale.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' multiDBDimPlot_vary_cells("CD14", object = "pbmc", cells.use.meta = "ident")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' multiDBDimPlot_vary_cells("CD14", cells.use.meta = "ident")

multiDBDimPlot_vary_cells <- function(var, object = DEFAULT,
                                      cells.use.meta,
                                      cells.use.levels = meta.levels(cells.use.meta,object),
                                      all.cells.plot = T,
                                      show.legend = F,
                                      add.single.legend = T,
                                      ncol = 3, nrow = NULL,
                                      add.title=T, axes.labels=F,
                                      min = NULL, max = NULL,
                                      data.type = "normalized",
                                      OUT.List = F,
                                      main = "make",
                                      all.cells.main = "All Cells",
                                      ...){

  #Interpret axes.labels: If left as FALSE, set lab to NULL so they will be removed.
  # If set to TRUE, set it to "make".
  lab <- if(!axes.labels) {NULL} else {"make"}

  #Determine if var is continuous vs discrete
  continuous <- FALSE
  if ((is.gene(var[1], object)) | (is.numeric(var))){
    continuous <- TRUE
  }
  if (is.meta(var[1], object)){
    if (is.numeric(meta(var,object))){
      continuous <- TRUE
    }
  }

  #Set a range if the var represents continuous data, (is.gene() or typeof(meta)!=integer OR character).
  if(continuous){
    the.range <- range(var_OR_get_meta_or_gene(var,object, data.type))
    min <- ifelse(is.null(min), the.range[1], min)
    max <- ifelse(is.null(max), the.range[2], max)
  }

  #Case 1: If var and cells.use.meta are NOT equal:
  #Need to ensure that colors are consistent.
  if(var[1] != cells.use.meta){
    #If data is continuous, then setting the range above has done the job of ensuring color consistency.
    if (continuous){
      plots <- lapply(cells.use.levels, function(X) {
        DBDimPlot(var, object,
                  cells.use = meta(cells.use.meta,object) == X,
                  auto.title = add.title,
                  xlab = lab,
                  ylab = lab,
                  main = ifelse(main == "make", X, main),
                  min = min, max = max,
                  data.type = data.type,
                  ...) +
          theme(legend.position = ifelse(show.legend, "right", "none"))
      })
      names(plots) <- cells.use.levels
    } else { #If data is discrete, then we need to make sure the right colors are use in each plot!
      levels <- meta.levels(var,object)
      plots <- lapply(cells.use.levels, function(X) {
        in.this.plot <- (1:length(levels))[levels %in%
                                             levels(as.factor(meta(var,object)[meta(cells.use.meta,object)==X]))]
        DBDimPlot(var, object,
                  cells.use = meta(cells.use.meta,object) == X,
                  auto.title = add.title,
                  xlab = lab,
                  colors = in.this.plot,
                  ylab = lab,
                  main = ifelse(main == "make", X, main),
                  min = min, max = max,,
                  data.type = data.type,
                  ...) +
          theme(legend.position = ifelse(show.legend, "right", "none"))
      })
      names(plots) <- cells.use.levels
    }
  }
  #Case 2: If var and cells.use.meta are equal
  #Need to vary the color with each new plot.
  if(var[1] == cells.use.meta){
    levels <- meta.levels(cells.use.meta, object)
    plots <- lapply((1:length(levels))[levels %in% cells.use.levels], function(X) {
      DBDimPlot(var, object,
                cells.use = meta(cells.use.meta,object) == levels[X],
                auto.title = add.title,
                colors = X,
                xlab = lab,
                ylab = lab,
                main = ifelse(main == "make", levels[X], main),
                min = min, max = max,,
                data.type = data.type,
                ...) +
        theme(legend.position = ifelse(show.legend, "right", "none"))
    })
    names(plots) <- cells.use.levels
  }
  #Generate all.cells.plot and legend explanation
  all.plot <- DBDimPlot(var, object,
                        auto.title = add.title,
                        xlab = lab,
                        ylab = lab,
                        min = min, max = max,,
                        data.type = data.type,
                        main = all.cells.main,
                        ...)
  legend <- cowplot::ggdraw(cowplot::get_legend(all.plot))
  if (all.cells.plot){
    plots$all <- all.plot + theme(legend.position = ifelse(show.legend, "right", "none"))
  }
  if (add.single.legend){
    plots$legend <- legend
  }
  if(OUT.List){
    return(plots)
  } else {
    return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
  }
}

#################### Helper Functions ########################

#### is.meta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot.
#'
#' @param test               "potential.meta.data.name" in quotes. REQUIRED.
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns TRUE if there is a meta.data slot named 'test' (or for Seurat objects, if test = "ident", will also give TRUE because my meta() function knows how to handle meta("ident"))
#' @examples
#' pbmc <- Seurat::pbmc_small
#' is.meta("age", object = pbmc) #Returns FALSE because there is not a meta slot named "age"
#' is.meta("nUMI", object = pbmc) #Returns TRUE because there is a meta slot named "age"
#' get.metas(object = pbmc) #get.metas() will give a list of all the meta.data slot names.
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' is.meta("age")
#' get.metas()

is.meta <- function(test, object=DEFAULT){

  #Bypass for ident for Seurat objects
  if(test=="ident" & grepl("Seurat",classof(object))) {return(TRUE)}

  #For all other cases...
  test %in% get.metas(object)
}

#### is.gene: Is this the name of a gene in my dataset? ####
#' Tests if an input is the name of a gene recovered in the dataset.
#'
#' @param test               "potential.gene.name" in quotes. REQUIRED.
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns TRUE if there is a row in the objects' data slot named 'test'.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' is.gene("CD14", object = "pbmc") # TRUE
#' is.gene("CD4", pbmc) # FALSE
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' is.gene("CD14")

is.gene <- function(test, object=DEFAULT){

  if(typeof(object)=="character")
  {
    if (classof(object)=="Seurat.v2"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    }
    if (classof(object)=="Seurat.v3"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object)))))
    }
    if (classof(object)=="RNAseq"){
      return(test %in% rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    }
  } else {
    if (class(object)=="Seurat.v2"){
      return(test %in% rownames(object@raw.data))
    }
    if (classof(object)=="Seurat.v3"){
      return(test %in% rownames(object))
    }
    if (class(object)=="RNAseq"){
      return(test %in% rownames(object@counts))
    }
  }
}

#### get.metas: prints the names of all the metadata lists for the object ####
#' Returns the names of all meta.data slots in the object.
#'
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the names of all meta.data slots in the object.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' get.metas(object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' get.metas()

get.metas <- function(object=DEFAULT){

  if(typeof(object)=="character"){
    names(eval(expr = parse(text = paste0(object,"@meta.data"))))
  } else {names(object@meta.data)}
}

#### get.genes: prints the names of all the genes for a Seurat or RNAseq ####
#' Returns the names of all genes within the data slot of the object.
#'
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the names of all genes within the data slot of the object.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' get.genes(object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' get.genes()

get.genes <- function(object=DEFAULT){

  if (classof(object)=="Seurat.v2"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    } else {return(rownames(object@raw.data))}
  }
  if (classof(object)=="Seurat.v3"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object)))))
    } else {return(rownames(object@raw.data))}
  }
  if (classof(object)=="RNAseq"){
    if(typeof(object)=="character"){
      rownames(eval(expr = parse(text = paste0(object,"@counts"))))
    } else {rownames(object@counts)}
  }
}

#### meta: for extracting the values of a particular metadata for all cells/samples ####
#' Returns the values of a meta.data  for all cells/samples
#'
#' @param meta               quoted "meta.data" slot name = REQUIRED. the meta.data slot that should be retrieved.
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the values of a meta.data slot, or the ident (clustering) slot if "ident" was given and the object is a Seurat object.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' meta("res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta("res.1")

meta <- function(meta, object=DEFAULT){
  if(typeof(object)=="character"){
    if(meta=="ident"){
      if(packageVersion("Seurat") >= '3.0.0'){
        return(as.character(Idents(object = eval(expr = parse(text = paste0(object))))))
      } else {return(as.character(eval(expr = parse(text = paste0(object,"@ident")))))}
    }
    else{return(eval(expr = parse(text = paste0(object,"@meta.data$'",meta, "'"))))}
  } else {
    if(meta=="ident"){
      if(packageVersion("Seurat") >= '3.0.0'){
        return(as.character(Idents(object)))
      } else {return(as.character(object@ident))}
    }
    else{return(eval(expr = parse(text = paste0("object@meta.data$",meta))))}
  }
}

#### gene: for extracting the expression values of a particular gene for all cells/samples ####
#' Returns the values of a gene for all cells/samples
#'
#' @param gene               quoted "gene" name = REQUIRED. the gene whose expression data should be retrieved.
#' @param object             "name" of the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param data.type          Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), or "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation)? Default = "normalized"
#' @return Returns the values of a meta.data slot, or the ident (clustering) slot if "ident" was given and the object is a Seurat object.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' gene("CD14", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' gene("CD14")

gene <- function(gene, object=DEFAULT, data.type = "normalized"){

  #Set up data frame for establishing how to deal with RNAseq or Seurat-v2 objects
  target <- data.frame(RNAseq = c("@data","@counts","error_Do_not_use_scaled_for_RNAseq_objects", "@samples"),
                       Seurat.v2 = c("@data","@raw.data","@scale.data", "@cell.names"),
                       Seurat.v3 = c("nope", "counts", "scale.data", "nope"),
                       stringsAsFactors = F,
                       row.names = c("normalized","raw","scaled","sample.names"))

  #Distinct functions if data.type = "relative" or "normalized.to.max" compared to all other options:
  if(data.type == "relative" | data.type == "normalized.to.max"){
    if(data.type == "relative"){
      #For "relative", recursive call to grab the 'normalized'  data that then has scaling run on top of it.
      OUT <- as.numeric(scale(gene(gene, object, "normalized")))
    } else {
      #For "normalized.to.max", recursive call to grab the 'normalized' data, then divide by its max.
      OUT <- gene(gene, object, "normalized")/max(gene(gene, object, "normalized"))
    }
    #For all other data.type options...
  } else {
    if (classof(object)!="Seurat.v3"){
      OUT <- eval(expr = parse(text = paste0(object,
                                             target[data.type,classof(object)],
                                             "[gene,",
                                             object,
                                             target["sample.names",classof(object)],
                                             "]")))
      #Change from sparse form if sparse
      OUT <- as.numeric(OUT)
      #Add names
      names(OUT) <- eval(expr = parse(text = paste0(object, target["sample.names",classof(object)])))
    } else {
      #Go from "object" to the actual object if given in character form
      object <- eval(expr = parse(text = paste0(object)))
      #Obtain expression
      if(data.type == "normalized"){
        OUT <- Seurat::GetAssayData(object)[gene,]
      } else {
        OUT <- Seurat::GetAssayData(object, slot = target[data.type,classof(object)])[gene,]
      }
      #Change from sparse form if sparse
      OUT <- as.numeric(OUT)
      #Add cellnames
      names(OUT) <- colnames(object)
    }
  }
  OUT
}

#### meta.levels: for obtaining the different classifications of a meta.data
#' Gives the distinct values of a meta.data slot (or ident)
#'
#' @param meta               quoted "meta.data.slot" name = REQUIRED. the meta.data slot whose potential values should be retrieved.
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param table.out          TRUE/FALSE. Default = FALSE. Whether the numbers of incidences of each level are wanted in addition to the level names themselves.
#' @return Returns the distinct values of a meta.data slot, or ident (clustering) slot if "ident" was given and the object is a Seurat object.  Can also return the counts of each as well.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' meta.levels("res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta.levels("res.1")

meta.levels <- function(meta, object = DEFAULT, table.out = F){

  if (table.out){
    table(meta(meta, object))
  } else {
    levels(as.factor(meta(meta, object)))
  }
}

#' Retrieves the proper slot of data
#'
#' @param data.type          "raw", "normalized", or "scaled". REQUIRED. which type of data is requested
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Given "raw", "normalized", or "scaled", this function will output the proper slot of a seurat or RNAseq object.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' which_data("normalized", "pbmc")

which_data <- function(data.type, object=DEFAULT){
  #Set up data frame for establishing how to deal with different input object types
  target <- data.frame(RNAseq = c("@data","@counts","error_Do_not_use_scaled_for_RNAseq_objects"),
                       Seurat.v2 = c("@data","@raw.data","@scale.data"),
                       Seurat.v3 = c("nope", "counts", "scale.data"),
                       stringsAsFactors = F,
                       row.names = c("normalized","raw","scaled"))
  if(classof(object)!="Seurat.v3"){
    #For RNAseq or Seurat-v2
    eval(expr = parse(text = paste0(object,
                                    target[data.type,classof(object)]
                                    )))
  } else {
    #For Seurat-v3
    #Go from "object" to the actual object if given in character form
    object <- eval(expr = parse(text = paste0(object)))
    #Obtain expression
    if(data.type == "normalized"){
      OUT <- Seurat::GetAssayData(object)
    } else {
      OUT <- Seurat::GetAssayData(object, slot = target[data.type,classof(object)])
    }
  }
}

#' Retrieves all the cells/samples
#'
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Given a seurat or RNAseq object, will return the cell.names or samples slot.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' all_cells("pbmc")

all_cells <- function(object = DEFAULT){
  object <- S4_2string(object)
  target <- data.frame(use = c("@samples", "@cell.names"),
                       row.names = c("RNAseq", "Seurat.v2"))
  if (classof(object)=="Seurat.v3"){
    return(colnames(x = eval(expr = parse(text = paste0(object)))))
  } else {
    eval(expr = parse(text = paste0(object, target[classof(object),])))
  }
}

#' Retrieves the list of cells.names to use
#'
#' @param cells.use either a logical or a list of names.
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Given a logical or a list of names (or NULL) will output the list of cells names.  For retrieval / standardization.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' which_cells(meta("ident","pbmc")=="0", "pbmc")

which_cells <- function(cells.use, object = DEFAULT){
  all.cells <- all_cells(object)
  if (is.null(cells.use)){
    return(all.cells)
  }
  if (is.logical(cells.use)){
    OUT <- all.cells[cells.use]
  } else {
    OUT <- cells.use
  }
  OUT
}


#' Outputs a string of gene expression / meta.data information for ploty hover display
#'
#' @param data.hover The data needed in the text output. = A list of metadata names, genes, or "ident", in quotes.
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param data.type For when grabbing gene expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @return Given a list of data to grab in data.hover, outputs the 'data name': data, 'data name': data, ... for every cell of the object
#' @examples
#' pbmc <- Seurat::pbmc_small
#' make_hover_strings(c("CD34","ident","non-genes/metas-will-be-ignored"), "pbmc", "normalized")

make_hover_strings <- function(data.hover, object, data.type = "normalized"){
  #Overall: if do.hover=T and data.hover has a list of genes / metas,
  # then for all cells, make a string "var1: var1-value\nvar2: var2-value..."
  if (is.null(data.hover)) {
    hover.string <- "Add gene or metadata \n names to hover.data"
  } else {
    features.info <- data.frame(row.names = all_cells(object))
    fill <- sapply(seq_len(length(data.hover)), function(i)
      (is.meta(data.hover[i],object) | is.gene(data.hover[i],object) | (data.hover[i]=="ident")))
    features.info <- sapply(seq_len(length(data.hover))[fill], function(i)
      features.info[,dim(features.info)[2]+1] <-
        var_OR_get_meta_or_gene(data.hover[i],object, data.type))
    names(features.info) <- data.hover[fill]
    hover.string <- sapply(seq_len(nrow(features.info)), function(row){
      paste(as.character(sapply(seq_len(ncol(features.info)), function(col){
        paste0(names(features.info)[col],": ",features.info[row,col])})),collapse = "\n")
    })
  }
  hover.string
}

#' Turns an S4 object into it's name in string form
#'
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return mainly for standardization within DittoSeq functions, outputs the string name and ensures objects can be handled in that form.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' S4_2string(pbmc)
#' S4_2string("pbmc")

S4_2string <- function(object = DEFAULT){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){
    object <- deparse(substitute(object))
  }
  object
}

#' Outputs gene expression data, metadata data, or clustering data
#'
#' @param var name of a metadata, gene, or "ident". = the data that should be grabbed
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param data.type For when extracting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @return determines what type of var is given, and outputs the gene expression data, metadata data, or clustering data.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' var_OR_get_meta_or_gene("CD14", "pbmc", "normalized")
#' var_OR_get_meta_or_gene("ident", "pbmc")
#' var_OR_get_meta_or_gene("nUMI", "pbmc")

var_OR_get_meta_or_gene <- function(var, object = DEFAULT, data.type){
  OUT <- var
  if(length(var)==1 & typeof(var)=="character"){
    #If "ident" pull the @ident object from the seurat object
    if(var == "ident"){OUT <- meta(var, object)}
    #If "is.meta" pull the @meta.data$"var" from the RNAseq or seurat object
    if(is.meta(var, object)){OUT <- meta(var, object)}
    #If "is.gene" pull the gene expression data from the RNAseq or seurat object
    if(is.gene(var, object)){OUT <- gene(var, object, data.type)}
    #Otherwise, var is likely a full set of data already, so just make Y = var
  }
  names(OUT) <- all_cells(object)
  OUT
}

#### extDim: for extracting the PC/tsne/IC/dim.reduction loadings of all cells/samples ####
#' Extracts the PC/tsne/IC/dim.reduction loadings of each cell/sample
#'
#' @param reduction.use      quoted "reduction" name. = REQUIRED. Common types are "pca", "tsne", "ica", "cca", "cca.aligned".
#' @param dim                #. which component to extract the loadings for (PC1 vs 2 vs 3)
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns a list where [[1]]=embeddings = the loadings of cell/sample for a given dimensional reduction component, and [[2]] = name = the string name that should be used to refer to the components if going into a plot.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' out <- extDim("pca", 1, object = "pbmc")
#' out$embeddings #will be the cell loadings
#' out$name # for a "pca" reduction will typically be "PC#", which is the proper way of
#'          # calling principal component 1, PC1, in a plot axis.
#' #If you wanted to plot the PC1 loadings of all cells in the different clusters...
#'
#' DBPlot(var = extDim("pca", 1, object = "pbmc")$embeddings,
#'        object = "pbmc",
#'        group.by = "res.1",
#'        color.by = "res.1",
#'        ylab = extDim("pca", 1, object = "pbmc")$name)
#'

extDim <- function(reduction.use, dim=1, object=DEFAULT){

  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){
    object <- deparse(substitute(object))
  }

  # If object is a Seurat object
  if (grepl("Seurat",classof(object))){
    if ("dr" %in% slotNames(eval(expr = parse(text = paste0(object))))){
      OUT <- list(eval(expr = parse(text = paste0(object,"@dr$",reduction.use,"@cell.embeddings[,",dim,"]"))))
      OUT[2] <- paste0(eval(expr = parse(text = paste0(object,"@dr$",reduction.use,"@key"))),dim)
    } else {
      OUT <- list(eval(expr = parse(text = paste0(object,"@reductions$",reduction.use,"@cell.embeddings[,",dim,"]"))))
      OUT[2] <- paste0(eval(expr = parse(text = paste0(object,"@reductions$",reduction.use,"@key"))),dim)
    }
  }

  if (classof(object)=="RNAseq"){
    OUT <- list(eval(expr = parse(text = paste0(object,"@reductions$",reduction.use,"$x[,",dim,"]"))))
    OUT[2] <- paste0(gen.key(reduction.use),dim)
  }
  names(OUT) <- c("embeddings","name")
  OUT
}

#### gen.key: for generating the proper axes label for a dimensional reduction type.
#' Creates a 'key' for properly refering to a dimensional reduction component when the full name (ex: PCA1) woul dbe wrong.
#'
#' @param reduction.use      quoted "reduction" name. = REQUIRED. Common types are "pca", "tsne", "ica", "cca", "cca.aligned".
#' @return Returns a string name that should be used to refer to the components if going into a plot.
#' @examples
#' gen.key("pca_5000_genes") # The function sees that pca is in the reduction name, so it outputs "PC"
#' #Output: "PC"
#' gen.key("tsne_5000_genes") # The function sees that pca is in the reduction name, so it outputs "PC"
#' #Output: "tSNE_"
#'

gen.key <- function (reduction.use){
  key <- reduction.use
  if (grepl("pca", reduction.use)){key <- "PC"}
  if (grepl("cca", reduction.use)){key <- "CC"}
  if (grepl("cca.aligned", reduction.use)){key <- "aligned.CC"}
  if (grepl("ica", reduction.use)){key <- "IC"}
  if (grepl("tsne", reduction.use)){key <- "tSNE_"}
  key
}

#### classof: for determining if 'object' is a Seurat or RNAseq ####
#' Returns the class of an object when given the name of the object in "quotes"
#'
#' @param object      quoted "object" name
#' @return Returns a the string name of the object's type.
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#' gen.key("pbmc")
#' #Output: "seurat"

classof <- function (object = DEFAULT){
  #DittoSeq works with the object in "string" form to work with it's DEFAULT setting method.

  #if object in "string" form, convert to the actual object
  if (typeof(object)=="character"){
    object <- eval(expr = parse(text = paste0(object)))
  }

  #class <- class(object)
  class <- class(object)

  #If a Seurat, need to add what version
  if (grepl("Seurat|seurat",class)){
    if(object@version >= '3.0.0') {
      #Then needs to be
      class <- "Seurat.v3"
    } else {
      class <- "Seurat.v2"
    }
  }
  class
}

###### grab_legend: Extract the legend from a ggplot ############
#' Extract the legend from a ggplot
#'
#' @description This function will extract the legend of any ggplot.
#' @param ggplot The ggplot object that you would like to grab the legend from.  Can be any of the single plots from DittoSeq except DBHeatmap and demux.calls.summary
#' @return The legend of a ggplot plot
#' @examples
#' #Grab data
#' pbmc <- Seurat::pbmc_small
#'
#' #Make a plot
#' DBDimPlot(pbmc)
#'
#' #Extract the legend:
#' grab_legend(DBDimPlot(pbmc))
#'
#' #Extract the legend of a stored plot:
#' ggplot <- DBDimPlot(pbmc)
#' grab_legend(ggplot)
grab_legend <- function(ggplot){
  cowplot::ggdraw(cowplot::get_legend(ggplot))
}

###### remove_legend: Remove the legend from a ggplot ############
#' Remove the legend from a ggplot
#'
#' @description This function will remove the legend of any ggplot.
#' @param ggplot The ggplot object that you would like to eliminate the legend from.  Can be any of the single plots from DittoSeq except DBHeatmap and demux.calls.summary
#' @return A ggplot plot with its legend removed.
#' @examples
#' #Grab data
#' pbmc <- Seurat::pbmc_small
#'
#' #Make a plot
#' DBDimPlot(pbmc)
#'
#' #Remove the legend:
#' remove_legend(DBDimPlot(pbmc))
#'
#' #Remove the legend of a stored plot:
#' ggplot <- DBDimPlot(pbmc)
#' remove_legend(ggplot)
remove_legend <- function(ggplot){
  ggplot + theme(legend.position = "none")
}

############################################################################################################

#### COLOR PANEL MODIFICATIONS BLOCK ####

############################################################################################################

#### Darken: For darkening colors ####
#' Darkens input colors by a set amount
#'
#' @description A wrapper for the darken function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, MYcolors.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = T.
#' @return Return a darkened version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Darken("blue") #"blue" = "#0000FF"
#' #Output: "#0000BF"
#' Darken(MYcolors[1:8]) #Works for multiple color inputs as well.

Darken <- function(colors, percent.change = 0.25, relative = T){

  colorspace::darken(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

#### Lighten: For lightening colors ####
#' Lightens input colors by a set amount
#'
#' @description A wrapper for the lighten function of the colorspace package.
#' @param colors the color(s) input. Can be a list of colors, for example, MYcolors.
#' @param percent.change # between 0 and 1. the percentage to darken by. Defaults to 0.25 if not given.
#' @param relative TRUE/FALSE. Whether the percentage should be a relative change versus an absolute one. Default = T.
#' @return Return a lighter version of the color in hexadecimal color form (="#RRGGBB" in base 16)
#' @examples
#' Lighten("blue") #"blue" = "#0000FF"
#' #Output: "#4040FF"
#' Lighten(MYcolors[1:8]) #Works for multiple color inputs as well.

Lighten <- function(colors, percent.change = 0.25, relative = T){

  colorspace::lighten(colors, amount = percent.change, space = "HLS", fixup = TRUE, method = ifelse(relative,"relative","absolute"))
}

#### Simulate: For simulating what a plot would look like as seen by a colorblind person ####
#' Simulates what a colorblind person would see for any DittoSeq plot!
#'
#' @description Essentially a wrapper function for colorspace's deutan(), protan(), and tritan() functions. This function will output any DittoSeq plot as it might look to an individual with one of the common forms of colorblindness: deutanopia/deutanomaly, the most common, is when the cones mainly responsible for red vision are defective. Protanopia/protanomaly is when the cones mainly responsible for green vision are defective. In tritanopia/tritanomaly, the defective cones are responsible for blue vision. Note: there are more severe color deficiencies that are even more rare. Unfortunately, for these types of color vision deficiency, only non-color methods, like lettering or shapes, will do much to help.
#' @param type The type of colorblindness that you want to simulate for. Options: "deutan", "protan", "tritan". Anything else, and you will get an error.
#' @param plot.function The plotting function that you want to use/simulate. not quoted. and make sure to remove the () that R will try to add.
#' @param color.panel The set of colors to be used.  Not required to be given, as contrary to the look of this documentation, it will still default to MYcolors when not provided.
#' @param ... other paramters that can be given to DittoSeq plotting functions, including color.panel, used in exactly the same way they are used for those functions. (contrary to the look of this documentation, color.panel will still default to MYcolors when not provided.)
#' @return Outputs any DittoSeq plot as it might look to a colorblind individual.
#' @examples
#' pbmc <- Seurat::pbmc_small
#' Simulate("deutan", DBDimPlot, var = "res.1", object = "pbmc", size = 2)
#' Simulate("protan", DBDimPlot, var = "res.1", object = "pbmc", size = 2)
#' Simulate("tritan", DBDimPlot, var = "res.1", object = "pbmc", size = 2)

Simulate <- function(type = "deutan", plot.function, color.panel = NULL, ...){

  #Check that type was given properly
  if(!(type=="deutan"|type=="protan"|type=="tritan")){
    return("type must be 'deutan', 'protan', or 'tritan'")
  }
  #Set the color panel
  if(is.null(color.panel)){
    color.p <- eval(expr = parse(text = paste0("colorspace::",type,"(MYcolors)")))
  } else {
    color.p <- eval(expr = parse(text = paste0("colorspace::",type,"(",color.panel,")")))
  }
  #Make the plot!
  plot.function(color.panel = color.p, ... )
}

################################################################################################################

#### Bulk RNA-Seq block ####

################################################################################################################

# The code in this block creates an object for DESeq2 data with similarish structure to a Seurat object.
# = helpful for analyzing bulk and single cell data side-by-side
# = helpful for analyzing bulk data if you are used to Seurat data structure
# = adds compatibility to my plotting and helper functions for bulk RNAseq data analyzed with DESeq2

#' The RNAseq Class
#'
#' @description The RNAseq object stores data analyzed in DESeq2 in a structure similar to a Seurat-object.  This is the data structure required for DittoSeq plottign functions to access bulk RNAseq data.  All that is needed to create an RNAseq object is a DESeqDataSet output from the DESeq() function.
#' @slot counts a matrix. The raw genes x samples counts data. It is recommended, but not required, that one of these should be given when a new RNBAseq object is created.
#' @slot dds a DESeqDataSet. The output of having run DESeq() on your data.
#' @slot data a matrix. The regularized log correction of the counts data generated by a call to DESeq's rlog function.
#' @slot meta.data a data.frame that contains meta-information about each sample. Autopopulated from the DESeq upon import, or added to manually afterward. Can be sample names, conditions, timepoints, Nreads (the number of reads).
#' @slot reductions a list of dimensional reductions that have been run. Any type of reduction can technically be supported, but the only one included with my package is pca by prcomp. Embeddings are stored in the reductions$pca$x slot.
#' @slot var.genes vector of genes with high coefficient of variation that passed an expression filter, and would therefore be included the default principal components analysis calculation.
#' @slot samples a vector of names of the samples.
#' @slot exp.filter a logical vector showing whether each gene passed the expression filter (default: at least 1 count in 75 percent of samples from each condition)
#' @slot CVs a numeric vector showing the coefficient of variation (mean divided by sd) for each gene in the dataset.
#' @slot other A great place to store any other data associated with your bulk experiment that does not fit elsewhere in the object.  Components of the list can be of any type.  Left empty by default, and is not altered or used by any of the DittoSeq functions.

Class <- setClass("RNAseq",
                  representation(
                    counts = "matrix",
                    dds = "ANY",
                    data = "matrix",
                    meta.data = "data.frame",
                    reductions = "list",
                    var.genes = "character",
                    samples = "character",
                    exp.filter = "logical",
                    CVs = "numeric",
                    other = "list"
                  ))

#### import.DESeq2 builds an RNAseq object with a DESeq input.  Can run PCA as well. ####
#' Creates an RNAseq object from a DESeq object.
#'
#' @description The first step of visualization of DESeq-analyzed bulk RNAseq data is running this function. Doing will extract meta.data information from the DESeq object, data using the rlog function, and if run_PCA=T, will populate all other slots as well.
#' @param dds                The output of running DESeq() on your data. = The DESeq2 object for your data. REQUIRED.
#' @param counts             Matrix. The raw counts data matrix.  Not required but HIGHLY RECOMMENDED.
#' @param run_PCA            TRUE/FALSE. Default is False. If set to true, prcomp PCA calculation will be carried out with the PCAcalc function.  var.genes, reductions$pca, exp.filter, and CVs slots will then all be populated. For more info, run ?PCAcalc
#' @param pc.genes           NULL or vector of genes. Alternately to the method of genes selection used by PCAcalc by default, a set of genes can be given here.  Default = NULL. If left that way, a per condition expression filter will be applied, followed by a selection of Ngenes number of genes that have the highest coefficient of variation (CV=mean over sd).
#' @param Ngenes             #. How many genes to use for the PCA calculation. (This number will ultimately be the length of the var.genes slot)
#' @param blind              TRUE/FALSE. Whether rlog estimation should be blinded to sample info. Run `?rlog` for more info about whether it should.  Defaults to TRUE, but that is NOT the correct way for all experiments.
#' @param percent.samples    # between 0 and 1. The percent of samples within each condition that must express the gene in order for a gene to be included in the PCA calculation.
#' @return Outputs an RNAseq object.
#' @examples
#'
#' #Generate mock RNAseq counts and a DESeq object from the mock data
#' # count tables from RNA-Seq data
#' counts.table <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#' colnames(counts.table) <- paste0("Sample",1:10)
#' conditions <- factor(rep(1:2, each=5))
#' # object construction
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(counts.table, DataFrame(conditions), ~ conditions)
#' dds <- DESeq(dds)
#'
#' # Recommended usage
#' # obj <- import.DESeq2(dds, counts = counts.table, run_PCA = TRUE)
#' #NOTE: the PCA calculation fails on this fake data because of it was normally randomized.
#' # Minimal input:
#' obj <- import.DESeq2(dds)

import.DESeq2 <- function(dds, #A DESeq object, *the output of DESeq()*
                          run_PCA = FALSE,#If changed to TRUE, function will:
                          # auto-populate var.genes, CVs, and pca fields.
                          pc.genes = NULL,
                          Ngenes = 2500, #How many genes to use for running PCA, (and how many genes
                          # will be stored in @var.genes)
                          blind = FALSE, #Whether or not the rlog estimation should be blinded to sample info.
                          # Run `?rlog` for more info
                          counts = NULL, #Raw Counts data, matrix with columns = genes and rows = samples.
                          # not required, but can be provided.
                          percent.samples = 75
){

  ########## Create the Object ########################
  #Create the object with whatever inputs were given, a.k.a. creates objects@counts and any other level
  # within str(object).  Will all be NULL unless provided in the function call
  object <- new("RNAseq", dds = dds)

  ########## Run Autopopulations ######################
  # Will run by default because this function requires a dds object to be given.
  ##populate dds
  object@dds <- dds
  ##populate @counts
  #Use the data provided if it was, otherwise, grab from the dds
  if (!(is.null(counts))) {object@counts <- counts
  } else {object@counts <- counts(dds)}
  ##populate @samples
  object@samples <- colnames(object@counts)
  ##populate some of @meta.data
  #   1st add samples, then Nreads.
  object@meta.data <- data.frame(Samples = object@samples,
                                 Nreads = colSums(object@counts))
  rownames(object@meta.data) <- object@samples

  ##Also add colData from dds to @meta.data slot
  #Turn colData into a data.frame, and merge that with current meta.data, BUT do not include any
  # dublicate sets.  For example, Samples will be ignored in colData because it was already grabbed
  # from the counts matrix
  object@meta.data <- cbind(object@meta.data,
                            data.frame(object@dds@colData@listData)[(!duplicated(
                              c(names(data.frame(object@dds@colData@listData)),names(object@meta.data)),
                              fromLast=T
                            ))[1:length(object@dds@colData@listData)]])

  ##populate data
  object@data <- SummarizedExperiment::assay(DESeq2::rlog(object@dds, blind = blind))

  ########## Will run if run_PCA = TRUE ##################
  ##Will populate: pca, CVs, and var.genes
  if (run_PCA){
    object <- PCAcalc(object = object,
                      genes.use = pc.genes,
                      Ngenes = Ngenes,
                      percent.samples = percent.samples,
                      name = "pca")
  }
  #OUTPUT: (This is how functions "work" in R.  The final line is what they return.)
  object
}

###### PCAcalc: For running the PCA calculation on an RNAseq object ############
#' Wrapper for running prcomp on an RNAseq object
#'
#' @description This function will run prcomp PCA calculation on either a set of genes given by the genes.use slot, or on the Ngenes that have the highest coefficient of variation (CV=mean/sd) after a per experimental condition (extracted from the dds) expression filter is applied.
#' @param object             the RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param genes.use          NULL or a vector of genes.  This will set the genes to be used by the prcomp PCA calculation.  NOTE: this list will not be used to populate the var.genes slot, but you can do that manually if you want to.
#' @param Ngenes             #. How many genes to use for the PCA calculation. (This number will ultimately be the length of the var.genes slot)
#' @param percent.samples    # between 0 and 100. The percent of samples within each condition that must express the gene in order for a gene to be included in the PCA calculation.
#' @param name               "name". Example: "pca". The name to be given to this pca reduction slot. -> redution$'name'
#' @return Outputs an RNAseq object with a new reductions$'name' slot.
#' @examples
#' # If bulkObject is the name of an RNAseq object in your workspace...
#' # PCAcalc("bulkObject", Ngenes = 2500, percent.samples = 75, name = "pca")
#' # Minimal input that does the same thing:
#' # PCAcalc("bulkObject")

PCAcalc <- function(object = DEFAULT,
                    genes.use = NULL,
                    Ngenes = 2500, #How many genes to use for running PCA, (and how many genes
                    # will be stored in @var.genes)
                    percent.samples = 75,
                    name = "pca"
){

  #Turn the percent.samples into a decimal named cutoff
  cutoff <- percent.samples/100

  #grab the object if given an object name
  if(typeof(object)=="character"){object <- eval(expr = parse(text = paste0(object)))}

  ######### IF no genes.use given, use the CVs and ExpFilter to pick genes ###########
  if(is.null(genes.use)){
    #Filter data to only the genes expressed in at least 75% of samples from test group (ONLY WORKS FOR ONE TEST GROUP)
    test_meta <- strsplit(as.character(object@dds@design), split = "~")[[2]][1]
    #Store this metadata as an easily accessible variable to speed up the next step.
    classes <- meta(test_meta, object)
    #For each gene, return TRUE if... the gene is expressed in >cutoff% of samples from each condition used to build the dds
    ##populate exp.filter
    object@exp.filter <- sapply(1:dim(object@counts)[1], function(X)
      #Ensures that the #of classifications = the number TRUEs in what the nested sapply produces
      length(levels(as.factor(classes)))==sum(
        #For each classification of the test variable, check for >= cutoff% expression accross samples
        #This half sets a variable to each of the saparate classifications,
        # and says to run the next lines on each
        sapply(levels(as.factor(classes)), function(Y)
          #This part of the function determine how many of the samples express the gene
          (sum(object@counts[X, classes==Y]>0))
          #This part compares the above to (the number of samples of the current classification*cutoff%)
          >= (sum(classes==Y)*cutoff)
        )
      )
    )
    data_for_prcomp <- as.data.frame(object@data)[object@exp.filter,]
    #calculate CV by dividing mean by sd
    ## populate CVs
    object@CVs <- apply(X = object@data, MARGIN = 1, FUN = sd)/apply(X = object@data, MARGIN = 1, FUN = mean)
    #Trim rlog data and RawCV_rlog by expression filter variable = object@exp.filter
    #arrange by CV_rank, higher CVs first
    data_for_prcomp<- data_for_prcomp[order(object@CVs[object@exp.filter], decreasing = T),]
    ##populate var.genes
    object@var.genes <- rownames(data_for_prcomp)[1:(min(Ngenes,dim(data_for_prcomp)[1]))]
    ##populate pca : Run PCA on the top Ngenes CV genes that survive the cutoff% expression per condition filter
    object@reductions$a1b2c3d57 <- prcomp(t(data_for_prcomp[1:Ngenes,]), center = T, scale = T)
  }
  ######### IF genes.use given, use them  ###########
  if(!(is.null(genes.use))){
    #Filter rlog to only the genes in genes.use
    data_for_prcomp <- as.data.frame(object@data)[genes.use,]
    ##populate pca : Run PCA on the top given genes.  DOES NOT USE THE Ngenes or percent.samples inputs!
    object@reductions$a1b2c3d57 <- list(prcomp(t(data_for_prcomp), center = T, scale = T))
  }

  names(object@reductions)[grep("a1b2c3d57",names(object@reductions))] <- name
  object
}
