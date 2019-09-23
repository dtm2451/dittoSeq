#### is.meta: Is this the name of a meta.data slot in my dataset? ####
#' Tests if an input is the name of a meta.data slot in a target object.
#'
#' @param test String or vector of strings, the "potential.metadata.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a meta.data slot within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed metadata slots if \code{return.values = TRUE}.
#' @details
#' For Seurat objects, also returns TRUE for the input \code{"ident"} because, for all dittoSeq visualiztions, \code{"ident"} will retrieve a Seurat objects' clustering slot.
#'
#' @seealso
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all metadata slots of an object
#' get.metas(pbmc)
#'
#' # To check if something is a metadata slot
#' is.meta("age", object = pbmc) # False
#' is.meta("nCount_RNA", object = pbmc) # True
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' is.meta("groups")
#'
#' # works for multiple test metas
#' is.meta(c("age","nCount_RNA"))
#'
#' # and with return.values = TRUE, returns all elements that are indeed slots
#' is.meta(c("age","nCount_RNA", "RNA_snn_res.0.8"),
#'     return.values = TRUE)
#'
#' @export
#' @import ggplot2
#' @importFrom utils packageVersion

is.meta <- function(test, object=DEFAULT, return.values=FALSE){
    if (typeof(object)=="S4") {
        object <- deparse(substitute(object))
    }
    if (return.values){
        return(test[is.meta(test, object, return.values=FALSE)])
    } else {
        metas <- get.metas(object)
        if (grepl("Seurat", .class_of(object))) {
            metas <- c(metas,"ident")
        }
        return(test %in% metas)
    }
}

#### is.gene: Is this the name of a gene in my dataset? ####
#' Tests if input is the name of a gene in a target object.
#'
#' @param test String or vector of strings, the "potential.gene.name"(s) to check for.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @param return.values Logical which sets whether the function returns a logical \code{TRUE}/\code{FALSE} versus the \code{TRUE} \code{test} values . Default = \code{FALSE}
#' REQUIRED, unless '\code{DEFAULT <- "object"}' has been run.
#' @return Returns a logical or logical vector indicating whether each instance in \code{test} is a gene within the \code{object}.
#' Alternatively, returns the values of \code{test} that were indeed genes if \code{return.values = TRUE}.
#' @seealso
#' \code{\link{get.genes}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' get.genes(pbmc)
#'
#' # To test if something is a gene in an object:
#' is.gene("CD14", object = "pbmc") # TRUE
#' is.gene("CD12345", pbmc) # FALSE
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped.
#' DEFAULT <- "pbmc"
#' is.gene("CD14")
#'   # TRUE
#'
#' # To test if many things are genes of an object
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"))
#'
#' # return.values input is especially useful in these cases.
#' is.gene(c("CD14", "IL32", "CD3E", "CD12345"), return.values = TRUE)
#'
#' @export

is.gene <- function(test, object=DEFAULT, return.values = FALSE){
    if (typeof(object)=="S4") {
        object <- deparse(substitute(object))
    }
    if (return.values) {
        return(test[is.gene(test, object, return.values=FALSE)])
    } else {
        return(test %in% get.genes(object))
    }
}

#### get.metas: prints the names of all the metadata lists for the object ####
#' Returns the names of all meta.data slots of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @return A string vector, returns the names of all metadata slots of the \code{object}
#' @seealso
#' \code{\link{is.meta}} for checking if certain metadata slots exist in an \code{object}
#'
#' \code{\link{meta}} for obtaining the contants of metadata slots
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all metadata slots of an object
#' get.metas(pbmc)
#'
#' @export

get.metas <- function(object=DEFAULT){
    if(.class_of(object)=="SingleCellExperiment"){
        if(typeof(object)=="character"){
            names(eval(expr = parse(text = paste0(object,"@colData"))))
        } else {names(object@colData)}
    } else {
        if(typeof(object)=="character"){
            names(eval(expr = parse(text = paste0(object,"@meta.data"))))
        } else {names(object@meta.data)}
    }
}

#### get.genes: prints the names of all the genes for a Seurat or RNAseq ####
#' Returns the names of all genes of a target object.
#'
#' @param object A target Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object, OR the name of the target object in "quotes".
#' @return A string vector, returns the names of all genes of the \code{object}.
#' @seealso
#' \code{\link{is.gene}} for returning all genes in an \code{object}
#'
#' \code{\link{gene}} for obtaining the expression data of genes
#'
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # To see all genes of an object
#' get.genes(pbmc)
#'
#' @export

get.genes <- function(object=DEFAULT){
    if (.class_of(object) %in% c("Seurat.v3","SingleCellExperiment")) {
        if (typeof(object)=="character") {
            return(rownames(eval(expr = parse(text = paste0(object)))))
        } else {
            return(rownames(object))
        }
    }
    if (.class_of(object)=="RNAseq") {
        if (typeof(object)=="character") {
            return(rownames(eval(expr = parse(text = paste0(object,"@counts")))))
        } else {
            return(rownames(object@counts))
        }
    }
    if (.class_of(object)=="Seurat.v2") {
        if (typeof(object)=="character") {
            return(rownames(eval(expr = parse(text = paste0(object,"@raw.data")))))
        } else {
            return(rownames(object@raw.data))
        }
    }
}

#### meta: for extracting the values of a particular metadata for all cells/samples ####
#' Returns the values of a meta.data  for all cells/samples
#'
#' @param meta String, the name of the "metadata" slot to grab. OR "ident" to retireve the clustering of a Seurat \code{object}.
#' @param object A Seurat, SingleCellExperiment, or \linkS4class{RNAseq} object to work with, OR the name of the object in "quotes".
#' @return Returns the values of a metadata slot, or the clustering slot if \code{meta = "ident"} and the \code{object} is a Seurat.
#' @seealso
#' \code{\link{meta.levels}} for returning just the unique discrete identities that exist within a metadata slot
#'
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{is.meta}} for testing whether something is the name of a metadata slot
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#' meta("RNA_snn_res.1", object = "pbmc")
#'
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta("RNA_snn_res.1")
#'
#' @export

meta <- function(meta, object=DEFAULT){
    if (typeof(object)=="S4") {
        object <- deparse(substitute(object))
    }

    if( meta=="ident" & grepl("Seurat", .class_of(object))) {
    # Retrieve clustering from Seurats
        if (grepl("v3",.class_of(object))) {
            return(as.character(Seurat::Idents(object =
                eval(expr = parse(text = paste0(
                    object))))))
        } else {
            return(as.character(
                eval(expr = parse(text = paste0(
                    object,"@ident")))))
        }
    }
    if (.class_of(object)=="SingleCellExperiment"){
        return(eval(expr = parse(text = paste0(
            object,"$'",meta,"'"))))
    } else {
        #RNAseq or Seurat non-ident
        return(eval(expr = parse(text = paste0(
            object,"@meta.data$'",meta, "'"))))
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
        OUT <- eval(expr = parse(text = paste0(
            "SingleCellExperiment::",
            target[data.type,.class_of(object)],
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
#' @param meta quoted "meta.data.slot" name = REQUIRED. the meta.data slot whose potential values should be retrieved.
#' @param object the Seurat or RNAseq object = REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param cells.use String vector of cells'/samples' names which should be included.
#' Alternatively, a Logical vector, the same length as the number of cells in the object, which sets which cells to include.
#' For the typically easier logical method, provide \code{USE} in \code{object@cell.names[USE]} OR \code{colnames(object)[USE]}).
#' @return Returns the distinct values of a metadata slot given to all cells/samples or for a subset of cells/samples.
#' (Alternatively, returns the distinct values of clustering if \code{meta = "ident"} and the object is a \code{Seurat} object).
#' @seealso
#' \code{\link{meta}} for returning an entire metadata slots of an \code{object}, not just the potential levels
#'
#' \code{\link{get.metas}} for returning all metadata slots of an \code{object}
#'
#' \code{\link{is.meta}} for testing whether something is the name of a metadata slot
#' @examples
#' library(Seurat)
#' pbmc <- Seurat::pbmc_small
#' meta.levels("RNA_snn_res.1", object = "pbmc")
#' # Note: if DEFAULT <- "pbmc" is run beforehand, the object input can be skipped completely.
#' DEFAULT <- "pbmc"
#' meta.levels("RNA_snn_res.1")
#' @export

meta.levels <- function(meta, object = DEFAULT, cells.use = NULL){
    if (typeof(object)=="S4"){
        object <- deparse(substitute(object))
    }
    meta.values <- as.character(meta(meta, object))
    if (!is.null(cells.use)){
        all.cells <- .all_cells(object)
        cells.use <- .which_cells(cells.use, object)
        meta.values <- meta.values[all.cells %in% cells.use]
    }
    levels(as.factor(meta.values))
}
