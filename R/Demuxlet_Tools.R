#' Imports 10X and demuxlet data and creates a Seurat object
#'
#' @param CellRangerLocations list of locations where barcodes.tsv, genes.tsv, and matrix.mtx are located.
#' @param Lane.names how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param Demux.best the location of the .best output file from running of demuxlet
#' @param verbose whether to print messages about the stage of this process that is currently being run.
#' @return Imports CellRanger 10X data from multiple lanes, merges those into a single Seurat object, then extracts the linked sample calls and demuxlet statistics from the demuxlet .best output file.  All relevant demuxlet information is stored as meta.data slots.
#' @examples
#' #Data required for an example would be rather large.  For an example, see the online vignette.

Import10XDemux <- function(CellRangerLocations, Lane.names=NULL, Demuxlet.best, verbose = TRUE){
  #CellRangerLocations = list of locations where barcodes.tsv, genes.tsv, and matrix.mtx are located.
  #Lane.names = how the lanes should be named, if you want to give them something different from Lane1, Lane2, Lane3...
  #Demux.best = the .best output file from running of demuxlet

  #If CellRangerLocations has length greater than one, create individual Seurat's for each and iteratively bring in each set.
  if(verbose){print("Building Seurat Object...")}
  if(length(CellRangerLocations)==1){
    OUT <- Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[1]))
  } else {
    #Create Lane.names if not given
    if(is.null(Lane.names)){Lane.names <- paste0("Lane",seq_len(length(CellRangerLocations)))}
    #Make Seurat Object
    if(verbose){print(paste0("  Making and merging ",Lane.names[1], " and ", Lane.names[2]))}
    OUT <- MergeSeurat(object1 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[1])),
                       add.cell.id1 = Lane.names[1],
                       object2 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[2])),
                       add.cell.id2 = Lane.names[2],
                       do.normalize = F)
    if(length(CellRangerLocations)>=3){
      for (X in 3:length(CellRangerLocations)){
        if(verbose){print(paste0("  Adding ", Lane.names[X]))}
        OUT <- MergeSeurat(object1 = OUT,
                           object2 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[X])),
                           add.cell.id2 = Lane.names[X],
                           do.normalize = F)
      }
    }
    if(verbose){print("  Adding 'Lane' information as meta.data")}
    OUT@meta.data$Lane <- sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[1]])
  }

  if(verbose){print("Extracting the Demuxlet information")}
    DEMUX.raw <- read.csv(Demuxlet.best, header=TRUE, sep="\t")
    DEMUX.call <- data.frame(doublet=sapply(as.character(DEMUX.raw$BEST),
                                            function(X) strsplit(X,'-')[[1]][[1]]),
                             sample=sapply(as.character(DEMUX.raw$BEST),
                                           function(x){strsplit(x,"-")[[1]][[2]]}),
                             barcode=sapply(as.character(DEMUX.raw$BARCODE),
                                            function(x){strsplit(x,"-")[[1]][[1]]}))
  if(verbose){print("  Matching barcodes")}
    # Strip barcodes from the cell names of the Seurat object
    cells <- sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[length(strsplit(X, "_")[[1]])]])
    # Remove any uneccesary cells from the Demux output (= cells not in this Seurat object)
    trim <- DEMUX.call[DEMUX.call$barcode %in% cells,]
    # Determine the indices to match Seurat cells to Demux matrix.
    inds <- sapply(cells,function(X) match(X,trim$barcode))
  if(verbose){print("  Adding Demuxlet as meta.data")}
    OUT@meta.data$Sample <- as.character(trim$sample[inds])
    OUT@meta.data$demux.doublet.call <- trim$doublet[inds]
    trim.out <- DEMUX.raw[DEMUX.call$barcode %in% cells,]
    OUT@meta.data$demux.RD.TOTL <- trim.out$RD.TOTL[inds]
    OUT@meta.data$demux.RD.PASS <- trim.out$RD.PASS[inds]
    OUT@meta.data$demux.RD.UNIQ <- trim.out$RD.UNIQ[inds]
    OUT@meta.data$demux.N.SNP <- trim.out$N.SNP[inds]
    OUT@meta.data$demux.PRB.DBL <- trim.out$PRB.DBL[inds]
  if(verbose){print("  Checking for barcode duplicates")}
    # Check if there are barcode duplicates accross 10X lanes.
    demux.barcode.dup <- array(F, dim = length(OUT@cell.names))
    demux.barcode.dup[duplicated(inds, fromLast=F)|duplicated(inds, fromLast=T)] <- T
    if (sum(demux.barcode.dup)>0){
      print("Warning: Cell barcodes were duplicated accross lanes of this dataset and may have been therefore artificially called as doublets by demuxlet. Recommended: Modify how your demuxlet was run to account for this.")
      OUT@meta.data$demux.barcode.dup <- demux.barcode.dup
    }
  OUT
}

# Test <- ImportMulti10XDemux(CellRangerLocations = c("/Users/danielbunis/Box Sync/Layering Analysis/10X/Raw-cellranger/CD4_matrix/",
#                                             "/Users/danielbunis/Box Sync/Layering Analysis/10X/Raw-cellranger/CD4-8_matrix/",
#                                             "/Users/danielbunis/Box Sync/Layering Analysis/10X/Raw-cellranger/CD8_matrix/",
#                                             "/Users/danielbunis/Box Sync/Layering Analysis/10X/Raw-cellranger/HSPC2_matrix/"),
#                     Lane.names = c("CD4","CD4.8","CD8","HSPC2"),
#                     Demuxlet.best = "../../pre-R/Demux/demux11.best")
# DEFAULT <- "Test"
# get.metas()
# meta.levels("Sample")
# meta.levels("Sample", table.out = T)
# meta.levels("Lane", table.out = T)
# table(meta("Lane"),meta("Sample"))


