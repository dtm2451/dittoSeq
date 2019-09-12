#### is.meta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot.
#'
#' @param test               "potential.meta.data.name" in quotes. REQUIRED.
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns TRUE if there is a meta.data slot named 'test' (or for Seurat objects, if test = "ident", will also give TRUE because my meta() function knows how to handle meta("ident"))
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' is.meta("age", object = pbmc) #Returns FALSE because there is not a meta slot named "age"
#' is.meta("nCount_RNA", object = pbmc) #Returns TRUE because there is a meta slot named "nCount_RNA"
#' get.metas(object = pbmc) #get.metas() will give a list of all the meta.data slot names.
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' is.meta("age")
#' get.metas()
#' @export
#' @import ggplot2
#' @importFrom utils packageVersion

is.meta <- function(test, object=DEFAULT){

  #Bypass for ident for Seurat objects
  if(test=="ident" & grepl("Seurat",.class_of(object))) {return(TRUE)}

  #For all other cases...
  test %in% get.metas(object)
}

#### is.gene: Is this the name of a gene in my dataset? ####
#' Tests if an input is the name of a gene recovered in the dataset.
#'
#' @param test "potential.gene.name" in quotes. REQUIRED.
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param value TRUE/FALSE. Default = FALSE.whether to return the gene names instead of TRUE/FALSE.
#' @return Returns TRUE if there is a row in the objects' data slot named 'test'.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#'
#' #For testing an individual query:
#' is.gene("CD14", object = "pbmc")
#'   # TRUE
#' is.gene("CD12345", pbmc)
#'   # FALSE, CD12345 is not a human gene.
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped.
#' DEFAULT <- "pbmc"
#' is.gene("CD14")
#'   # TRUE
#'
#' #The function works for sets of gene-queries as well
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"))
#'   #TRUE TRUE TRUE FALSE
#' # value input is especially useful in these cases.
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"), value = TRUE)
#'  #"CD14" "IL32" "CD3E"
#' @export

is.gene <- function(test, object=DEFAULT, value = FALSE){
  if (value){
    #Return names of included genes
    return(test[is.gene(test, object, value=FALSE)])
  } else {
    #Return TRUE/FALSE
    return(test %in% get.genes(object))
  }
}

#### get.metas: prints the names of all the metadata lists for the object ####
#' Returns the names of all meta.data slots in the object.
#'
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the names of all meta.data slots in the object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' get.metas(object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' get.metas()
#' @export

get.metas <- function(object=DEFAULT){
  if(.class_of(object)=="SingleCellExperiment"){
    #SingleCellExperiment
    if(typeof(object)=="character"){
      names(eval(expr = parse(text = paste0(object,"@colData"))))
    } else {names(object@colData)}
  } else {
    #Non- SingleCellExperiment
    if(typeof(object)=="character"){
      names(eval(expr = parse(text = paste0(object,"@meta.data"))))
    } else {names(object@meta.data)}
  }
}

#### get.genes: prints the names of all the genes for a Seurat or RNAseq ####
#' Returns the names of all genes within the data slot of the object.
#'
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the names of all genes within the data slot of the object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' get.genes(object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' get.genes()
#' @export

get.genes <- function(object=DEFAULT){
  if (.class_of(object)=="Seurat.v2"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    } else {return(rownames(object@raw.data))}
  }
  if (.class_of(object)=="Seurat.v3"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object)))))
    } else {return(rownames(object))}
  }
  if (.class_of(object)=="RNAseq"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    } else {return(rownames(object@counts))}
  }
  if (.class_of(object)=="SingleCellExperiment"){
    if(typeof(object)=="character"){
      return(return(rownames(counts(eval(expr = parse(text = paste0(object)))))))
    } else {return(rownames(counts(object)))}
  }
}

#### get.reductions: prints the names of all the dimensional reductions that have been run ####
#' Returns the names of dimensional reduction slots within the object.
#'
#' @param object             the SingleCellExperiment, Seurat, or RNAseq object = REQUIRED unless `DEFAULT <- "object"` has been run.
#' @return Returns the names of all dimensional reduction slots within the object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' get.reductions(object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' get.reductions()
#' @export

get.reductions <- function(object=DEFAULT){
  if (.class_of(object)=="Seurat.v2"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@dr")))))
    } else {return(names(object@dr))}
  }
  if (.class_of(object)=="Seurat.v3"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@reductions")))))
    } else {return(names(object@reductions))}
  }
  if (.class_of(object)=="RNAseq"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@reductions")))))
    } else {return(names(object@reductions))}
  }
  if (.class_of(object)=="SingleCellExperiment"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@reducedDims")))))
    } else {return(names(object@reducedDims))}
  }
}

#### meta: for extracting the values of a particular metadata for all cells/samples ####
#' Returns the values of a meta.data  for all cells/samples
#'
#' @param meta               quoted "meta.data" slot name = REQUIRED. the meta.data slot that should be retrieved.
#' @param object             the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Returns the values of a meta.data slot, or the ident (clustering) slot if "ident" was given and the object is a Seurat object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' meta("RNA_snn_res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta("RNA_snn_res.1")
#' @export

meta <- function(meta, object=DEFAULT){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){
    object <- deparse(substitute(object))
  }
  #Name of object given in "quotes"
  if(meta=="ident" & grepl("Seurat", .class_of(object))){
    if(eval(expr = parse(text = paste0(object, "@version"))) >= '3.0.0'){
      #Seurat.v3, ident
      return(as.character(Seurat::Idents(object = eval(expr = parse(text = paste0(object))))))
    } else {
      #Seurat.v2, ident
      return(as.character(eval(expr = parse(text = paste0(object,"@ident")))))
    }
  } else {
    if (.class_of(object)=="SingleCellExperiment"){
      #SingleCellExperiment
      return(eval(expr = parse(text = paste0(object,"$'",meta,"'"))))
    } else {
      #RNAseq or Seurat non-ident
      return(eval(expr = parse(text = paste0(object,"@meta.data$'",meta, "'"))))
    }
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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' gene("CD14", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' gene("CD14")
#' @export

gene <- function(gene, object=DEFAULT, data.type = "normalized"){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4"){
    object <- deparse(substitute(object))
  }
  #Set up data frame for establishing how to deal with RNAseq or Seurat-v2 objects
  target <- data.frame(RNAseq = c("@data","@counts","error_Do_not_use_scaled_for_RNAseq_objects", "@samples"),
                       Seurat.v2 = c("@data","@raw.data","@scale.data", "@cell.names"),
                       Seurat.v3 = c("nope", "counts", "scale.data", "nope"),
                       SingleCellExperiment = c("logcounts", "counts", "error_do_not_use_scaled_for_SCE_objects", "nope"),
                       stringsAsFactors = FALSE,
                       row.names = c("normalized","raw","scaled","sample.names"))

  #Distinct functions if data.type = "relative", "normalized.to.max", or "raw.normalized.to.max" compared to all other options:
  if(data.type == "relative" | data.type == "normalized.to.max" | data.type == "raw.normalized.to.max"){
    if(data.type == "relative"){
      #For "relative", recursive call to grab the 'normalized' data, then has scale run on top of it.
      OUT <- as.numeric(scale(gene(gene, object, "normalized")))
    }
    if(data.type == "normalized.to.max"){
      #For "normalized.to.max", recursive call to grab the 'normalized' data, then divide by its max.
      OUT <- gene(gene, object, "normalized")/max(gene(gene, object, "normalized"))
    }
    if(data.type == "raw.normalized.to.max"){
      #For "raw.normalized.to.max", recursive call to grab the 'raw' data, then divide by its max.
      OUT <- gene(gene, object, "raw")/max(gene(gene, object, "raw"))
    }
    #For all other data.type options...
  } else {
    if (.class_of(object)!="Seurat.v3" & .class_of(object)!="SingleCellExperiment"){
      OUT <- eval(expr = parse(text = paste0(object,
                                             target[data.type,.class_of(object)],
                                             "[gene,",
                                             object,
                                             target["sample.names",.class_of(object)],
                                             "]")))
      #Change from sparse form if sparse
      OUT <- as.numeric(OUT)
      #Add names
      names(OUT) <- eval(expr = parse(text = paste0(object, target["sample.names",.class_of(object)])))
    } else {
      if (.class_of(object)=="Seurat.v3"){
        #Go from "object" to the actual object if given in character form
        object <- eval(expr = parse(text = paste0(object)))
        #Obtain expression
        if(data.type == "normalized"){
          OUT <- Seurat::GetAssayData(object)[gene,]
        } else {
          OUT <- Seurat::GetAssayData(object, slot = target[data.type,.class_of(object)])[gene,]
        }
      } else {
        #SingleCellExperiment
        OUT <- eval(expr = parse(text = paste0(target[data.type,.class_of(object)],
                                               "(", object, ")[gene,]")))
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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' meta.levels("RNA_snn_res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta.levels("RNA_snn_res.1")
#' @export

meta.levels <- function(meta, object = DEFAULT, table.out = FALSE){
  if (table.out){
    table(meta(meta, object))
  } else {
    levels(as.factor(meta(meta, object)))
  }
}

###### grab_legend: Extract the legend from a ggplot ############
#' Extract the legend from a ggplot
#'
#' @description This function will extract the legend of any ggplot.
#' @param ggplot The ggplot object that you would like to grab the legend from.  Can be any of the single plots from dittoSeq except DBHeatmap and demux.calls.summary
#' @return The legend of a ggplot plot
#' @examples
#' library(Seurat)
#'
#' #Grab data
#' pbmc <- Seurat::pbmc_small
#' DEFAULT <- "pbmc"
#'
#' #Make a plot
#' DBDimPlot("ident")
#'
#' #Extract the legend:
#' grab_legend(DBDimPlot("ident"))
#'
#' #Extract the legend of a stored plot:
#' ggplot <- DBDimPlot("ident")
#' grab_legend(ggplot)
#' @export
grab_legend <- function(ggplot){
  cowplot::ggdraw(cowplot::get_legend(ggplot))
}

###### remove_legend: Remove the legend from a ggplot ############
#' Remove the legend from a ggplot
#'
#' @description This function will remove the legend of any ggplot.
#' @param ggplot The ggplot object that you would like to eliminate the legend from.  Can be any of the single plots from dittoSeq except DBHeatmap and demux.calls.summary
#' @return A ggplot plot with its legend removed.
#' @examples
#' library(Seurat)
#'
#' #Grab data
#' pbmc <- Seurat::pbmc_small
#' DEFAULT <- "pbmc"
#'
#' #Make a plot
#' DBDimPlot("ident")
#'
#' #Remove the legend:
#' remove_legend(DBDimPlot("ident"))
#'
#' #Remove the legend of a stored plot:
#' ggplot <- DBDimPlot("ident")
#' remove_legend(ggplot)
#' @export
remove_legend <- function(ggplot){
  ggplot + theme(legend.position = "none")
}
