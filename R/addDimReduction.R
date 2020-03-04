#' Add a prcomp pca calculation to a SingleCellExperiment object containing bulk or single-cell data
#'
#' @param object the \code{\linkS4class{SingleCellExperiment}} object.
#' @param prcomp a prcomp output which will be added to the \code{object}
#' @param name String name for the reduction slot.
#' Normally, this will be "pca", but you can hold any number of PCA calculations so long as a unique \code{name} is given to each.
#' This will become the name of the slot and what should be provided to the \code{reduction.use} input when making a \code{\link{dittoDimPlot}}.
#' When the name given is the same as that of a slot that already exists inside the \code{object}, the previous slot is replaced with the newly provided data.
#' @param key String, like "PC", which sets the default axes-label prefix when this reduction is used for making a \code{\link{dittoDimPlot}}
#' @return Outputs an \code{\linkS4class{SingleCellExperiment}} object with an added or replaced pca reduction slot.
#' @seealso
#' \code{\link{addDimReduction}} for adding other types of dimensionality reductions
#'
#' \code{\link{importDittoBulk}} for initial import of bulk RNAseq data into dittoSeq as a \code{\linkS4class{SingleCellExperiment}}.
#'
#' \code{\link{dittoDimPlot}} for visualizing how samples group within added dimensionality reduction spaces
#' @examples
#'
#' example("importDittoBulk", echo = FALSE)
#'
#' # Calculate PCA with prcomp
#' #   NOTE: This is typically not done with all genes in a dataset.
#' #   The inclusion of this example code is not an endorsement of a particular
#' #   method of PCA. Consult yourself, a bioinformatician, or literature for
#' #   tips on proper techniques.
#' calc <- prcomp(t(logcounts(myRNA)), center = TRUE, scale = TRUE)
#'
#' myRNA <- addPrcomp(
#'     object = myRNA,
#'     prcomp = calc)
#'
#' # Now we can visualize conditions metadata on a PCA plot
#' dittoDimPlot(myRNA, "conditions", reduction.use = "pca", size = 3)
#' @author Daniel Bunis
#' @export

addPrcomp <- function(object, prcomp, name = "pca", key = "PC") {
    addDimReduction(object,prcomp$x,name,key)
}

#' Add any dimensionality reduction space to a SingleCellExperiment object containing bulk or single-cell data
#'
#' @param object the bulk or single-cell \code{\linkS4class{SingleCellExperiment}} object to add the dimensionality reduction to.
#' (dittoSeq utilizes the SingleCellExperiment object even for bulk data because it provides a convenient slots for all data that dittoSeq requires)
#' @param embeddings a numeric matrix or matrix-like object, with number of rows equal to ncol(object), containing the coordinates of all cells / samples within the dimensionality reduction space.
#' @param name String name for the reduction slot. Example: "pca".
#' This will become the name of the slot, and what should be provided to the \code{reduction.use} input when making a \code{\link{dittoDimPlot}}.
#' When the name given is the same as that of a slot that already exists inside the \code{object}, the previous slot is replaced with the newly provided data.
#' @param key String, like "PC", which sets the default axes-label prefix when this reduction is used for making a \code{\link{dittoDimPlot}}.
#' If nothing is provided, a key will be automatically generated.
#' @return Outputs a \code{\linkS4class{SingleCellExperiment}} object with an added or replaced dimensionality reduction slot.
#' @seealso
#' \code{\link{addPrcomp}} for a prcomp specific PCA import wrapper
#'
#' \code{\link{importDittoBulk}} for initial import of bulk RNAseq data into dittoSeq as a \code{\linkS4class{SingleCellExperiment}}.
#'
#' \code{\link{dittoDimPlot}} for visualizing how samples group within added dimensionality reduction spaces
#' @examples
#'
#' example("importDittoBulk", echo = FALSE)
#'
#' # Calculate PCA
#' #   NOTE: This is typically not done with all genes in the dataset.
#' #   The inclusion of this example code is not an endorsement of a particular
#' #   method of PCA. Consult yourself, a bioinformatician, or literature for
#' #   tips on proper techniques.
#' embeds <- prcomp(t(logcounts(myRNA)), center = TRUE, scale = TRUE)$x
#'
#' myRNA <- addDimReduction(
#'     object = myRNA,
#'     embeddings = embeds,
#'     name = "pca",
#'     key = "PC")
#'
#' # Visualize conditions metadata on a PCA plot
#' dittoDimPlot(myRNA, "conditions", reduction.use = "pca", size = 3)
#'
#' @author Daniel Bunis
#' @importFrom SingleCellExperiment reducedDim<-
#' @export

addDimReduction <- function(
    object, embeddings, name, key = .gen_key(name)) {

    # Add cellnames
    rownames(embeddings) <- colnames(object)
    # Add dim names based on key
    colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))
    # Add
    SingleCellExperiment::reducedDim(object, name) <- embeddings
    object
}
