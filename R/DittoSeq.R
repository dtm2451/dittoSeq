#### is.meta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot.
#' @importFrom utils packageVersion
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

is.meta <- function(test, object=DEFAULT){

  #Bypass for ident for Seurat objects
  if(test=="ident" & grepl("Seurat",classof(object))) {return(TRUE)}

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
  if(classof(object)=="SingleCellExperiment"){
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
  if (classof(object)=="Seurat.v2"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
    } else {return(rownames(object@raw.data))}
  }
  if (classof(object)=="Seurat.v3"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object)))))
    } else {return(rownames(object))}
  }
  if (classof(object)=="RNAseq"){
    if(typeof(object)=="character"){
      return(rownames(eval(expr = parse(text = paste0(object,"@counts")))))
    } else {return(rownames(object@counts))}
  }
  if (classof(object)=="SingleCellExperiment"){
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
  if (classof(object)=="Seurat.v2"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@dr")))))
    } else {return(names(object@dr))}
  }
  if (classof(object)=="Seurat.v3"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@reductions")))))
    } else {return(names(object@reductions))}
  }
  if (classof(object)=="RNAseq"){
    if(typeof(object)=="character"){
      return(names(eval(expr = parse(text = paste0(object,"@reductions")))))
    } else {return(names(object@reductions))}
  }
  if (classof(object)=="SingleCellExperiment"){
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
  if(typeof(object)=="character"){
    #Name of object given in "quotes"
    if(meta=="ident" & grepl("Seurat", classof(object))){
      if(eval(expr = parse(text = paste0(object, "@version"))) >= '3.0.0'){
        #Seurat.v3, ident
        return(as.character(Seurat::Idents(object = eval(expr = parse(text = paste0(object))))))
      } else {
        #Seurat.v2, ident
        return(as.character(eval(expr = parse(text = paste0(object,"@ident")))))
      }
    } else {
      if (classof(object)=="SingleCellExperiment"){
        #SingleCellExperiment
        return(eval(expr = parse(text = paste0(object,"$'",meta,"'"))))
      } else {
        #RNAseq or Seurat non-ident
        return(eval(expr = parse(text = paste0(object,"@meta.data$'",meta, "'"))))
      }
    }
  } else {
    # Actual object given
    if(meta=="ident"){
      if(packageVersion("Seurat") >= '3.0.0'){
        #Seurat.v3, ident
        return(as.character(Seurat::Idents(object)))
      } else {
        #Seurat.v2, ident
        return(as.character(object@ident))
      }
    } else {
      if (classof(object)=="SingleCellExperiment"){
        #SingleCellExperiment
        return(eval(expr = parse(text = paste0("object$'",meta,"'"))))
      } else {
        #RNAseq or Seurat non-ident
        return(eval(expr = parse(text = paste0("object@meta.data$",meta))))
      }
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

  #Set up data frame for establishing how to deal with RNAseq or Seurat-v2 objects
  target <- data.frame(RNAseq = c("@data","@counts","error_Do_not_use_scaled_for_RNAseq_objects", "@samples"),
                       Seurat.v2 = c("@data","@raw.data","@scale.data", "@cell.names"),
                       Seurat.v3 = c("nope", "counts", "scale.data", "nope"),
                       SingleCellExperiment = c("logcounts", "counts", "error_do_not_use_scaled_for_SCE_objects", "nope"),
                       stringsAsFactors = FALSE,
                       row.names = c("normalized","raw","scaled","sample.names"))

  #Distinct functions if data.type = "relative" or "normalized.to.max" compared to all other options:
  if(data.type == "relative" | data.type == "normalized.to.max" | data.type == "raw.normalized.to.max"){
    if(data.type == "relative"){
      #For "relative", recursive call to grab the 'normalized'  data that then has scaling run on top of it.
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
    if (classof(object)!="Seurat.v3" & classof(object)!="SingleCellExperiment"){
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
      if (classof(object)=="Seurat.v3"){
        #Go from "object" to the actual object if given in character form
        object <- eval(expr = parse(text = paste0(object)))
        #Obtain expression
        if(data.type == "normalized"){
          OUT <- Seurat::GetAssayData(object)[gene,]
        } else {
          OUT <- Seurat::GetAssayData(object, slot = target[data.type,classof(object)])[gene,]
        }
      } else {
        #SingleCellExperiment
        OUT <- eval(expr = parse(text = paste0(target[data.type,classof(object)],
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

#' Retrieves the proper slot of data
#'
#' @param data.type          "raw", "normalized", or "scaled". REQUIRED. which type of data is requested
#' @param object             the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Given "raw", "normalized", or "scaled", this function will output the proper slot of a seurat or RNAseq object.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' which_data("normalized", "pbmc")

which_data <- function(data.type, object=DEFAULT){
  #Set up data frame for establishing how to deal with different input object types
  target <- data.frame(RNAseq = c("@data","@counts","error_Do_not_use_scaled_for_RNAseq_objects"),
                       Seurat.v2 = c("@data","@raw.data","@scale.data"),
                       Seurat.v3 = c("nope", "counts", "scale.data"),
                       SingleCellExperiment = c("logcounts", "counts", "error_do_not_use_scaled_for_SCE_objects"),
                       stringsAsFactors = FALSE,
                       row.names = c("normalized","raw","scaled"))
  if(classof(object)=="SingleCellExperiment"){
    OUT <- as.matrix(eval(expr = parse(text = paste0(target[data.type,classof(object)],
                                                     "(", object, ")" ))))
  } else {
    if(classof(object)!="Seurat.v3"){
      #For RNAseq or Seurat-v2
      OUT <- eval(expr = parse(text = paste0(object,
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
  OUT
}

#' Retrieves all the cells/samples
#'
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @return Given a seurat or RNAseq object, will return the cell.names or samples slot.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' all_cells("pbmc")

all_cells <- function(object = DEFAULT){
  object <- S4_2string(object)
  target <- data.frame(use = c("@samples", "@cell.names"),
                       row.names = c("RNAseq", "Seurat.v2"))
  if (classof(object)=="Seurat.v3" | classof(object)=="SingleCellExperiment"){
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
#' library(Seurat)
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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' make_hover_strings(c("CD34","ident","non-genes/metas-will-be-ignored"), "pbmc", "normalized")

make_hover_strings <- function(data.hover, object, data.type = "normalized"){
  #Overall: if do.hover=TRUE and data.hover has a list of genes / metas,
  # then for all cells, make a string "var1: var1-value\nvar2: var2-value..."
  if (is.null(data.hover)) {
    hover.string <- "Add gene or metadata \n names to hover.data"
  } else {
    features.info <- data.frame(row.names = all_cells(object))
    fill <- sapply(seq_along(data.hover), function(i)
      (is.meta(data.hover[i],object) | is.gene(data.hover[i],object) | (data.hover[i]=="ident")))
    features.info <- sapply(seq_along(data.hover)[fill], function(i)
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
#' library(Seurat)
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

#' Outputs the given data, gene expression data, metadata data, or clustering data
#'
#' @param var name of a metadata, gene, or "ident". = the data that should be grabbed
#' @param object the Seurat or RNAseq object to draw from = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param data.type For when extracting expression data: Should the data be "normalized" (data slot), "raw" (raw.data or counts slot), "scaled" (the scale.data slot of Seurat objects), "relative" (= pulls normalized data, then uses the scale() function to produce a relative-to-mean representation), or "normalized.to.max" (= pulls normalized data, then divides by the maximum value)? DEFAULT = "normalized"
#' @return determines what type of var is given, and outputs var itself or the gene expression data, metadata data, or clustering data refered to.
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' var_OR_get_meta_or_gene("CD14", "pbmc", "normalized")
#' var_OR_get_meta_or_gene("ident", "pbmc")
#' var_OR_get_meta_or_gene("nCount_RNA", "pbmc")
#' @export

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
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' out <- extDim("pca", 1, object = "pbmc")
#' out$embeddings #will be the cell loadings
#' out$name # for a "pca" reduction will typically be "PC#", which is the proper way of
#'          # calling principal component 1, PC1, in a plot axis.
#' #If you wanted to plot the PC1 loadings of all cells in the different clusters...
#'
#' DBPlot(var = extDim("pca", 1, object = "pbmc")$embeddings,
#'        object = "pbmc",
#'        group.by = "RNA_snn_res.1",
#'        color.by = "RNA_snn_res.1",
#'        ylab = extDim("pca", 1, object = "pbmc")$name)
#'
#' @export

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

  if (classof(object)=="SingleCellExperiment"){
    OUT <- list(eval(expr = parse(text = paste0(object,"@reducedDims$",reduction.use,"[,",dim,"]"))))
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
  if (grepl("pca|PCA", reduction.use)){key <- "PC"}
  if (grepl("cca|CCA", reduction.use)){key <- "CC"}
  if (grepl("cca.aligned", reduction.use)){key <- "aligned.CC"}
  if (grepl("ica|ICA", reduction.use)){key <- "IC"}
  if (grepl("tsne|tSNE|TSNE", reduction.use)){key <- "tSNE_"}
  key
}

#### classof: for determining if 'object' is a Seurat or RNAseq ####
#' Returns the class of an object when given the name of the object in "quotes"
#'
#' @param object quoted "object" name
#' @return Returns a the string name of the object's type.
#' @examples
#' library(Seurat)
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
#' @param ggplot The ggplot object that you would like to eliminate the legend from.  Can be any of the single plots from DittoSeq except DBHeatmap and demux.calls.summary
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
