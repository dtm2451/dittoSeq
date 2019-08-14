################################# dittoSingleAxisDataGather ####################################
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
    extra.vars = NULL, cells.use = NULL, data.type = "normalized", do.hover = FALSE, data.hover = c(var, extra.vars)){
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
  names <- c(main.var, "grouping", "color")
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
