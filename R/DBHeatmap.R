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
#' @param highlight.genes = A list of genes whose names you would like to show, c("Gene1","Gene2","Gene3").  Only these genes will be named in the resulting heatmap.
#' @param ... other arguments that will be passed to pheatmap. For scRNAseq data, I recommend running with cluster_cols=FALSE first because clustering of thousands of cells can tak a long time!  Common ones: main for adding a title, cluster_cols or cluster_rows = FALSE to turn clustering off for cells or genes, treeheight_row or treeheight_col = 0 to remove or a positive # for setting how large the trees on the side/top should be drawn
#' @description This function is a wrapper for the pheatmap function of the pheatmap package.  Given a set of genes, it will spit out a basic heatmap from your RNAseq data, or just certain cells/samples.
#' @return Makes a heatmap where
#' library(Seurat)
#' pbmc <- pbmc_small
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
                      data.out=FALSE, highlight.genes = NULL,
                      ...){

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

  #Make a labels_row input for displaying only certain genes if genes were given to highlight.genes input
  Row_names <- NULL
  if(!(is.null(highlight.genes))){
    if(sum(highlight.genes %in% genes)>0){
      #Subset to highlight.genes in the genelist
      highlight.genes <- highlight.genes[highlight.genes %in% genes]
      #Make a Row_names variable with the rownames of the data.
      Row_names <- rownames(data)
      #Overwrite all non-highlight genes rownames to ""
      Row_names[-(match(highlight.genes,Row_names))] <- ""
    }
  }

  #Establish cell/sample names
  Col_names <- NULL
  if(!(is.null(cell.names.meta))){
    Col_names <- as.character(meta(cell.names.meta, object)[all.cells %in% cells.use])
  }

  if(data.out){
    OUT <- list(data,heatmap.colors,Col_annot,Col_annot_colors, Col_names,Row_names,
                paste0("pheatmap(mat = data, color = heatmap.colors, annotation_col = Col_annot, ",
                       "annotation_colors = Col_annot_colors, scale = 'row', ",
                       "labels_col = Col_names, labels_row = Row_names"))
    names(OUT)<- c("data","heatmap.colors","Col_annot","Col_annot_colors","Col_names","Row_names","pheatmap.script")
    return(OUT)
  }

  #Make the Heatmap
  HM <- pheatmap::pheatmap(mat = data,
                           color = heatmap.colors,
                           annotation_col = Col_annot,
                           annotation_colors = Col_annot_colors,
                           scale = "row",
                           labels_col = Col_names,
                           labels_row = Row_names,
                           ...)

  #DONE: Return the heatmap
  return(HM)
}
