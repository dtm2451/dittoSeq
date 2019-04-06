#' Imports 10X and demuxlet data and creates a Seurat object
#'
#' @param CellRangerLocations list of locations where barcodes.tsv, genes.tsv, and matrix.mtx are located.
#' @param Lane.names how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param Demuxlet.best the location of the .best output file from running of demuxlet
#' @param verbose whether to print messages about the stage of this process that is currently being run.
#' @return Imports CellRanger 10X data from multiple lanes, merges those into a single Seurat object, then extracts the linked sample calls and demuxlet statistics from the demuxlet .best output file.  All relevant demuxlet information is stored as meta.data slots.
#' @examples
#' #Data required for an example would be rather large.  For an example, see the online vignette.

Import10XDemux <- function(CellRangerLocations, Lane.names=NA, Demuxlet.best, verbose = TRUE){
  #CellRangerLocations = list of locations where barcodes.tsv, genes.tsv, and matrix.mtx are located.
  #Lane.names = how the lanes should be named, if you want to give them something different from Lane1, Lane2, Lane3...
  #Demux.best = the .best output file from running of demuxlet

  #Create Lane.names if not given
  if(is.na(Lane.names[1])){Lane.names <- paste0("Lane",seq_len(length(Seurats)))}

  #If CellRangerLocations has length greater than one, create individual Seurat's for each and iteratively bring in each set.
  if(verbose){print("Building Seurat Object...",quote=F)}
  if(length(CellRangerLocations)==1){
    OUT <- Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[1]))
  } else {
    #Make Seurat Object
    if(verbose){print(paste0("  Making and merging ",Lane.names[1], " and ", Lane.names[2]),quote=F)}
    OUT <- Seurat::MergeSeurat(object1 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[1])),
                       add.cell.id1 = Lane.names[1],
                       object2 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[2])),
                       add.cell.id2 = Lane.names[2],
                       do.normalize = F)
    if(length(CellRangerLocations)>=3){
      for (X in 3:length(CellRangerLocations)){
        if(verbose){print(paste0("  Adding ", Lane.names[X]),quote=F)}
        OUT <- Seurat::MergeSeurat(object1 = OUT,
                           object2 = Seurat::CreateSeuratObject(Seurat::Read10X(CellRangerLocations[X])),
                           add.cell.id2 = Lane.names[X],
                           do.normalize = F)
      }
    }
    if(verbose){print("  Adding 'Lane' information as meta.data",quote=F)}
    OUT@meta.data$Lane <-
      if(length(CellRangerLocations)>1){
        sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[1]])
      } else {
        Lane.names[1]
      }
  }

  if(verbose){print("Extracting the Demuxlet information",quote=F)}
    DEMUX.raw <- read.csv(Demuxlet.best, header=TRUE, sep="\t")
    DEMUX.call <- data.frame(doublet=sapply(as.character(DEMUX.raw$BEST),
                                            function(X) strsplit(X,'-')[[1]][[1]]),
                             sample=sapply(as.character(DEMUX.raw$BEST),
                                           function(x){strsplit(x,"-")[[1]][[2]]}),
                             barcode=sapply(as.character(DEMUX.raw$BARCODE),
                                            function(x){strsplit(x,"-")[[1]][[1]]}))
  if(verbose){print("  Matching barcodes",quote=F)}
    # Strip barcodes from the cell names of the Seurat object
    cells <- sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[length(strsplit(X, "_")[[1]])]])
    # Remove any uneccesary cells from the Demux output (= cells not in this Seurat object)
    trim <- DEMUX.call[DEMUX.call$barcode %in% cells,]
    # Determine the indices to match Seurat cells to Demux matrix.
    inds <- sapply(cells,function(X) match(X,trim$barcode))
  if(verbose){print("  Adding Demuxlet as meta.data",quote=F)}
    OUT@meta.data$Sample <- as.character(trim$sample[inds])
    OUT@meta.data$demux.doublet.call <- trim$doublet[inds]
    trim.out <- DEMUX.raw[DEMUX.call$barcode %in% cells,]
    OUT@meta.data$demux.RD.TOTL <- trim.out$RD.TOTL[inds]
    OUT@meta.data$demux.RD.PASS <- trim.out$RD.PASS[inds]
    OUT@meta.data$demux.RD.UNIQ <- trim.out$RD.UNIQ[inds]
    OUT@meta.data$demux.N.SNP <- trim.out$N.SNP[inds]
    OUT@meta.data$demux.PRB.DBL <- trim.out$PRB.DBL[inds]
  if(verbose){
    print("  Checking for barcode duplicates",quote=F)
    print("",quote=F)
    print("",quote=F)}
  # Check if there are barcode duplicates accross 10X lanes.
  demux.barcode.dup <- array(F, dim = length(OUT@cell.names))
  demux.barcode.dup[duplicated(inds, fromLast=F)|duplicated(inds, fromLast=T)] <- T
  if (sum(demux.barcode.dup)>0){
    print("Warning: Cell barcodes were duplicated accross lanes of this dataset and may have been therefore artificially called as doublets by demuxlet. Recommended: Modify how your demuxlet was run to account for this.",quote=F)
    OUT@meta.data$demux.barcode.dup <- demux.barcode.dup
    print("",quote=F)
    print("",quote=F)
  }
  print("SUMMARY:",quote=F)
  print(paste0(length(CellRangerLocations),
               " CellRanger outputs were combined into a Seurat object with demuxlet data extracted into metadata."),quote=F)
  print(paste0("The average number of SNPs per cell was: ",round(mean(OUT@meta.data$demux.N.SNP),1)),quote=F)
  print(paste0("Out of ",length(OUT@cell.names)," cells total, Demuxlet assigned:"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="SNG", na.rm = TRUE), " cells or ", round(100*sum(OUT@meta.data$demux.doublet.call=="SNG")/length(OUT@cell.names),1), "% as singlets"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="DBL", na.rm = TRUE), " cells or ", round(100*sum(OUT@meta.data$demux.doublet.call=="DBL")/length(OUT@cell.names),1), "% as doublets"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="AMB", na.rm = TRUE), " cells as too ambiguous to call."),quote=F)
  OUT
}

#' Imports demuxlet data into a (merged) Seurat object
#'
#' @param Seurats The set of Seurat objects that represent your dataset (each is assumed to represent individual 10X lanes).
#' @param Lane.names how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param Demuxlet.best the location of the .best output file from running of demuxlet
#' @param verbose whether to print messages about the stage of this process that is currently being run.
#' @return Imports Seurat objects, each expected to represent a single lane of droplet generation (or separate pool of samples for another sc-seq method, or sequencing lane, etc.), merges those into a single Seurat object, then extracts the linked sample calls and demuxlet statistics from the demuxlet .best output file.  All relevant demuxlet information is stored as meta.data slots.
#' @examples
#' #Data required for an example would be rather large.  For an example, see the online vignette.

ImportSeuratDemux <- function(Seurats, Lane.names=NA, Demuxlet.best, verbose = TRUE){
  #Seurats = list of Seurat objects (each representing individual 10X lanes).
  #Lane.names = how the lanes should be named, if you want to give them something different from Lane1, Lane2, Lane3...
  #Demux.best = the .best output file from running of demuxlet

  #Create Lane.names if not given
  if(is.na(Lane.names[1])){Lane.names <- paste0("Lane",seq_len(length(Seurats)))}

  #If Seurats has length greater than one, merge individual Seurat's iteratively.
  if(verbose){print("Merging Seurat Objects...",quote=F)}
  if(length(Seurats)==1){
    OUT <- Seurats
  } else {
    #Make Seurat Object
    if(verbose){print(paste0("  Making and merging ",Lane.names[1], " and ", Lane.names[2]),quote=F)}
    OUT <- Seurat::MergeSeurat(object1 = Seurats[[1]],
                               add.cell.id1 = Lane.names[1],
                               object2 = Seurats[[2]],
                               add.cell.id2 = Lane.names[2],
                               do.normalize = F)
    if(length(Seurats)>=3){
      for (X in 3:length(Seurats)){
        if(verbose){print(paste0("  Adding ", Lane.names[X]),quote=F)}
        OUT <- Seurat::MergeSeurat(object1 = OUT,
                                   object2 = Seurats[[X]],
                                   add.cell.id2 = Lane.names[X],
                                   do.normalize = F)
      }
    }
    if(verbose){print("  Adding 'Lane' information as meta.data",quote=F)}
    OUT@meta.data$Lane <-
      if(length(Seurats)>1){
        sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[1]])
      } else {
        Lane.names[1]
      }
  }

  if(verbose){print("Extracting the Demuxlet information",quote=F)}
  DEMUX.raw <- read.csv(Demuxlet.best, header=TRUE, sep="\t")
  DEMUX.call <- data.frame(doublet=sapply(as.character(DEMUX.raw$BEST),
                                          function(X) strsplit(X,'-')[[1]][[1]]),
                           sample=sapply(as.character(DEMUX.raw$BEST),
                                         function(x){strsplit(x,"-")[[1]][[2]]}),
                           barcode=sapply(as.character(DEMUX.raw$BARCODE),
                                          function(x){strsplit(x,"-")[[1]][[1]]}))
  if(verbose){print("  Matching barcodes",quote=F)}
  # Strip barcodes from the cell names of the Seurat object
  cells <- sapply(OUT@cell.names, function(X) strsplit(X, "_")[[1]][[length(strsplit(X, "_")[[1]])]])
  # Remove any uneccesary cells from the Demux output (= cells not in this Seurat object)
  trim <- DEMUX.call[DEMUX.call$barcode %in% cells,]
  # Determine the indices to match Seurat cells to Demux matrix.
  inds <- sapply(cells,function(X) match(X,trim$barcode))
  if(verbose){print("  Adding Demuxlet as meta.data",quote=F)}
  OUT@meta.data$Sample <- as.character(trim$sample[inds])
  OUT@meta.data$demux.doublet.call <- trim$doublet[inds]
  trim.out <- DEMUX.raw[DEMUX.call$barcode %in% cells,]
  OUT@meta.data$demux.RD.TOTL <- trim.out$RD.TOTL[inds]
  OUT@meta.data$demux.RD.PASS <- trim.out$RD.PASS[inds]
  OUT@meta.data$demux.RD.UNIQ <- trim.out$RD.UNIQ[inds]
  OUT@meta.data$demux.N.SNP <- trim.out$N.SNP[inds]
  OUT@meta.data$demux.PRB.DBL <- trim.out$PRB.DBL[inds]
  if(verbose){
    print("  Checking for barcode duplicates",quote=F)
    print("",quote=F)
    print("",quote=F)}
  # Check if there are barcode duplicates accross 10X lanes.
  demux.barcode.dup <- array(F, dim = length(OUT@cell.names))
  demux.barcode.dup[duplicated(inds, fromLast=F)|duplicated(inds, fromLast=T)] <- T
  if (sum(demux.barcode.dup)>0){
    print("Warning: Cell barcodes were duplicated accross lanes of this dataset and may have been therefore artificially called as doublets by demuxlet. Recommended: Modify how your demuxlet was run to account for this.",quote=F)
    OUT@meta.data$demux.barcode.dup <- demux.barcode.dup
    print("",quote=F)
    print("",quote=F)
  }
  print("SUMMARY:",quote=F)
  print(paste0(length(Seurats),
               " Seurat objects were combined with demuxlet data extracted into metadata."),quote=F)
  print(paste0("The average number of SNPs per cell was: ",round(mean(OUT@meta.data$demux.N.SNP),1)),quote=F)
  print(paste0("Out of ",length(OUT@cell.names)," cells total, Demuxlet assigned:"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="SNG", na.rm = TRUE), " cells or ", round(100*sum(OUT@meta.data$demux.doublet.call=="SNG")/length(OUT@cell.names),1), "% as singlets"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="DBL", na.rm = TRUE), " cells or ", round(100*sum(OUT@meta.data$demux.doublet.call=="DBL")/length(OUT@cell.names),1), "% as doublets"),quote=F)
  print(paste0("     ",sum(OUT@meta.data$demux.doublet.call=="AMB", na.rm = TRUE), " cells as too ambiguous to call."),quote=F)
  OUT
}

#' Plots the number of SNPs sequenced per droplet
#'
#' @param object                 the Seurat Object = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param group.by               "metadata" to use for separating values. Default is "Lane".
#' @param color.by               "metadata" to use for coloring. Default is "Lane".
#' @param plots                  Default = c("jitter","boxplot"). = the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot". NOTE: The order matters, so use c("back","middle","front") when inputing multiple to put them in the order you want.
#' @param boxplot.color The color of the lines of the boxplot.
#' @param ...                    extra arguments passed to DBPlot
#' @return For a given Seurat object, plots a summary of how many SNPs were available to Demuxlet for making sample annotation calls.  Assumes that the number of SNPs are stored in a 'demux.N.SNPs' metadata slot, as would be the case if the Seurat object was created with the Import10XDemux function.
#' @examples
#' #Data required for an example would be rather large.  For an example, see the online vignette.
demux.SNP.summary <- function(object = DEFAULT, group.by = "Lane", color.by = "Lane", plots = c("jitter","boxplot"), boxplot.color = "grey30", ...){
  DBPlot("demux.N.SNP", object, group.by = group.by, color.by = color.by, plots = plots, boxplot.color = boxplot.color, ...)
}

#' Plots the number of annotations per sample, per lane
#'
#' @param object the Seurat Object = name of object in "quotes". REQUIRED, unless `DEFAULT <- "object"` has been run.
#' @param singlets.only Whether to only show data for cells called as singlets by demuxlet. Default is TRUE. Note: if doublets are included, only one of their sample calls will be used.
#' @param main plot title. Default = "Sample Annotations by Lane"
#' @param sub plot subtitle
#' @param ylab y axis label, default is "Annotations"
#' @param xlab x axis label, default is "Sample"
#' @param color bars color. Default is blue.
#' @param theme A complete ggplot theme. Default is a modified theme_bw().
#' @param rotate.labels whether sample names / x-axis labels should be rotated or not. Default is TRUE.
#' @return For a given Seurat object, summarizes how many cells in each lane were anotated to each sample.  Assumes that the Sample calls of each cells, and which lane each cell belonged to, are stored in 'Sample' and 'Lane' metadata slots, respectively, as would be the case if the Seurat object was created with the Import10XDemux function.
#' @examples
#' #Data required for an example would be rather large.  For an example, see the online vignette.
demux.calls.summary <- function(object = DEFAULT, singlets.only = TRUE,
                                main = "Sample Annotations by Lane", sub = NULL,
                                ylab = "Annotations", xlab = "Sample",
                                color = MYcolors[2], theme = NULL, rotate.labels = TRUE
){
  #Change object to character if not already
  object <- S4_2string(object)
  #Populate cells.use with a list of names, based on the singlets.only variable.
  cells.use <- if(singlets.only){
    which_cells(meta("demux.doublet.call",object)=="SNG", object)
  } else {
    all_cells(object)
  }
  #Establish the full list of cell/sample names
  all.cells <- all_cells(object)
  #Set theme
  if(is.null(theme)){
    theme <- theme_bw()+theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_blank(),
                              panel.border = element_blank())
  }

  dat <- as.data.frame.matrix(table(meta("Sample")[cells.use %in% all.cells],
                                    meta("Lane")[cells.use %in% all.cells]))
  dat$Sample <- row.names(dat)
  dat.m <- reshape2::melt(dat, "Sample")
  p <- ggplot(data = dat.m) +
    geom_col(aes(x = Sample, y = value), fill = color) +
    geom_hline(yintercept=0) +
    facet_wrap(~ variable, ncol = 1, scales = "fixed", strip.position = "left") +
    ggtitle(main, subtitle = sub) +
    xlab(xlab) + ylab(ylab) +
    theme
  #Rotate Labels if rotate.labels = TRUE
  if (rotate.labels) {p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1, size=12))}

  return(p)
}
