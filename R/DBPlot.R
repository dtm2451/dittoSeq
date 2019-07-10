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
#' @param main                   plot title. Default = "make", if left as make, a title will be automatically generated.  To remove, set to NULL.
#' @param sub                    plot subtitle
#' @param theme                  Allows setting of a theme. Default = theme_classic when nothing is provided.
#' @param ylab                   y axis label, defaults to "var" or "var expression" if var is a gene.  Set to NULL to remove.
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
#' @param boxplot.fill           whether the boxplot should be filled in or not.
#' @param reorder.x              sequence of numbers from 1:length(meta.levels(group.by)) for providing a new order for the samples.  Default = alphabetical then numerical.
#' @param legend.show            TRUE/FALSE. Whether the legend should be displayed. Default = TRUE.
#' @param title.legend           whether to leave the title for the plot's legend
#' @param vlnplot.lineweight sets the thickness of the line that outlines the violin plots.
#' @return Makes a plot where continuous data, grouped by sample, age, cluster, etc., on the x-axis is shown on the y-axis by a violin plot, boxplot, and/or dots (or other shapes)
#' @examples
#' library(Seurat)
#' pbmc <- pbmc_small
#' DBPlot("CD14", object = "pbmc", group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBPlot("CD14", group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")

DBPlot <- function(var, object = DEFAULT, group.by, color.by,
                   shape.by = "", cells.use = NULL, plots = c("jitter","vlnplot"),
                   data.type = "normalized", do.hover = FALSE, data.hover = var,
                   color.panel = MYcolors, colors = c(1:length(color.panel)),
                   theme = theme_classic(), main = "make", sub = NULL,
                   ylab = "make", y.breaks = NULL, min = NULL, max = NULL,
                   xlab = NULL, labels = NULL, rotate.labels = TRUE,
                   hline=NULL, hline.linetype = "dashed", hline.color = "black",
                   jitter.size=1, jitter.width=0.2, jitter.color = "black", jitter.shapes=c(16,15,17,23,25,8),
                   jitter.shape.legend.size = 3,
                   boxplot.width = 0.2, boxplot.color = "black", boxplot.show.outliers = NA, boxplot.fill =T,
                   vlnplot.lineweight = 1,
                   reorder.x = 1:length(meta.levels(group.by, object)),
                   legend.show = TRUE, title.legend = F){

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

  #Unless ylab has been changed from default ("make"), set default y-labels if var is "metadata" or "gene".
  # If set to "var" use the 'var'.  If set to 'NULL', ylab will be removed later.
  if (!(is.null(ylab))){
    if (ylab == "make" & length(var)==1 & is.character(var)) {
      if(is.meta(var, object)){ylab <- var}
      if(is.gene(var, object)){ylab <- paste0(var," expression")}
    }
    if(ylab=="var") {ylab <- var}
    if(ylab=="make") {ylab <- NULL}
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
      p <- p + geom_violin(size = vlnplot.lineweight)
    }
  }

  #Add horizontal lines if given.
  if (!is.null(hline))  {p <- p + geom_hline(yintercept=hline, linetype= hline.linetype, color = hline.color)}

  #Set default title if 'main' was unchanged from its default ("make") and is a gene or metadata.
  if (!(is.null(main))){
    if (main == "make" & is.character(var) & length(var)==1 & (is.meta(var, object) | is.gene(var, object))){
      main <- var
    }
    if (main == "make"){
      main <- NULL
    }
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

  #Remove legend, if warrented
  if (!legend.show) { p <- remove_legend(p) }

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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' genes <- c("CD8A","CD3E","FCER1A","CD14")
#' multiDBPlot(genes, object = "pbmc",
#'             group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object
#'       # input can be skipped completely.
#' DEFAULT <- "pbmc"
#' multiDBPlot(genes,
#'             group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' #To make it output a grid that is 2x2, to add y-axis labels
#' # instead of titles, and to show legends...
#' multiDBPlot(genes,
#'             group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
#'             nrow = 2, ncol = 2,              #Make it 2x2
#'             add.title = FALSE, ylab = TRUE,  #Add y axis labels instead of titles
#'             show.legend = TRUE)              #Show legends
#' # To eliminate the "expression", change ylab = TRUE to ylab = "var"
#' multiDBPlot(genes,
#'             group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
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

  plots <- lapply(vars, function(X) {
    DBPlot(X, object, group.by, color.by,
           ylab = if(typeof(ylab.input)=="character"){ifelse((ylab.input == "var"),X,return('Error: ylab must be TRUE/FALSE or "var"'))} else {if(ylab.input){"make"}else{NULL}},
           ...) + theme(legend.position = ifelse(show.legend, "right", "none"))
  })

  #Output
  if (OUT.List){
    names(plots) <- vars
    return(plots)
  } else {
    return(gridExtra::grid.arrange(grobs=plots, ncol = ncol, nrow = nrow))
  }
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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' genes <- get.genes("pbmc")[1:30]
#' DBPlot_multi_var_summary(genes, object = "pbmc",
#'                          group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object
#'       # input can be skipped completely.
#' DEFAULT <- "pbmc"
#' DBPlot_multi_var_summary(genes,
#'                          group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1")
#'
#' # To change it to have the violin plot in the back, a jitter on
#' #  top of that, and a white boxplot with no fill in front:
#' DBPlot_multi_var_summary(genes, object = "pbmc",
#'                          group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
#'                          plots = c("vlnplot","jitter","boxplot"),
#'                          boxplot.color = "white", boxplot.fill = FALSE)

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
      p <- p + ylim(min, max(dat$Y))
    }
    if (!(is.null(max))){
      p <- p + ylim(min(dat$Y), max)
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

  #Remove legend, if warrented
  if (!legend.show) { p <- remove_legend(p) }

  #DONE. Return the plot
  if(do.hover & ("jitter" %in% plots)){
    return(plotly::ggplotly(p, tooltip = "text"))
  } else {
    return(p)
  }
}
