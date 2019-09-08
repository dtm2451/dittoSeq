#' Extracts Demuxlet information into a pre-made Seurat object
#' @importFrom utils read.csv
#'
#' @param Seurat.name The "quoted" name of the Seurat object to add demuxlet information to.
#' @param Lane.info.meta A meta.data
#' @param Lane.names how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param Demuxlet.best the location of the .best output file from running of demuxlet
#' @param verbose whether to print messages about the stage of this process that is currently being run.
#' @param bypass.check for running this function anyway even when meta.data slots already around would be over-written.
#' @return Imports Demuxlet sample calls and statistics into a Seurat object.  All relevant demuxlet information is stored as meta.data slots.  Works in 2 differnet ways:  If a Lane.info.meta is given, this metadata is copied into a "Lane" metadata becuase that is where my summarizer functions look by default.  If no Lane.info.meta is given, information about lane origin is assumed to be in the form of "-#" at the end of the cellnames; the form they would be in if 10X lanes were combined with `cellranger aggr`.
#' @examples
#' #Data required for an example would be rather large.
#' # For an example, see the online vignette at:
#' # https://github.com/dtm2451/DittoSeq/tree/master/Demuxlet-Vignette
#' @export

ImportDemux2Seurat <- function(Seurat.name,
                               Lane.info.meta = NULL, Lane.names=NA,
                               Demuxlet.best,
                               bypass.check = FALSE,
                               verbose = TRUE){
  #Turn the object into a "name" if a full object was given
  if (typeof(Seurat.name)=="S4"){ Seurat.name <- deparse(substitute(Seurat.name)) }

  Seurat <- eval(expr = parse(text = Seurat.name))

  #Check if there are meta.data that would be over.written
  Check <- c("Lane", "Sample", "demux.doublet.call", "demux.RD.TOTL", "demux.RD.PASS",
             "demux.RD.UNIQ", "demux.N.SNP", "demux.PRB.DBL", "demux.barcode.dup")
  if (sum(Check %in% get.metas(Seurat.name))>0){
    if(bypass.check){
      print(paste0("Caution: '",
                   paste0(Check[Check %in% get.metas(Seurat.name)], collapse = "' and '"),
                   "' are going to bw overwritten."),
            quote = FALSE)
    } else {
      print("To proceed anyway, set `bypass.check = TRUE`.")
      print(paste0("WARNING: ",
                   paste0(Check[Check %in% get.metas(Seurat.name)], collapse = " and "),
                   " would be overwritten."),
            quote = FALSE)
      return(Seurat)
    }
  }
  #Create Lane.names if not given
  if(is.na(Lane.names[1])){
    Lane.names <- meta.levels(Lane.info.meta, Seurat.name)
  }
  # Obtain the cell.names
  cell.names <- .all_cells(Seurat.name)
  #Add lane metadata information
  if(verbose){print("Adding 'Lane' information as meta.data",quote=FALSE)}
  if(!(is.null(Lane.info.meta))){
    Seurat@meta.data$Lane <- meta(Lane.info.meta, Seurat.name)
    lane.idents <- as.factor(as.character(meta(Lane.info.meta, Seurat.name)))
  } else {
    # Make object with just the `#` of `-#` parts of cellnames
    lane.idents <- as.factor(c(sapply(cell.names, function(X){
      strsplit(X, split = "-")[[1]][length(strsplit(X, split = "-")[[1]])]
    })))
    # Create the Lane metadata
    Seurat@meta.data$Lane <- lane.idents
    # Change from # to Lane.names
    levels(Seurat@meta.data$Lane) <- Lane.names
  }
  Seurat@meta.data$Lane <- as.factor(as.character(Seurat@meta.data$Lane))

  #Extract Demuxlet data
  if(verbose){print("Extracting the Demuxlet information",quote=FALSE)}
  DEMUX.raw <- read.csv(Demuxlet.best, header=TRUE, sep="\t")
  DEMUX.call <- data.frame(doublet=sapply(as.character(DEMUX.raw$BEST),
                                          function(X) strsplit(X,'-')[[1]][[1]]),
                           sample=sapply(as.character(DEMUX.raw$BEST),
                                         function(x){strsplit(x,"-")[[1]][[2]]}),
                           barcode=as.character(DEMUX.raw$BARCODE))
  if(verbose){print("  Matching barcodes",quote=FALSE)}
  # Strip barcodes in cell.names from any of the extra info that may have been added by Seurat (normally "text_" at start of names)
  cells <- sapply(cell.names, function(X) strsplit(X, split = "_")[[1]][length(strsplit(X, split = "_")[[1]])])
  # Remove any uneccesary cells from the Demux output (= cells not in this Seurat object)
  trim <- DEMUX.call[DEMUX.call$barcode %in% cells,]
  # Determine the indices to match Seurat cells to Demux matrix.
  inds <- sapply(cells,function(X) match(X,trim$barcode))
  # Add demux info
  if(verbose){print("  Adding Demuxlet as meta.data",quote=FALSE)}
  Seurat@meta.data$Sample <- as.character(trim$sample[inds])
  Seurat@meta.data$demux.doublet.call <- trim$doublet[inds]
  trim.out <- DEMUX.raw[DEMUX.call$barcode %in% cells,]
  Seurat@meta.data$demux.RD.TOTL <- trim.out$RD.TOTL[inds]
  Seurat@meta.data$demux.RD.PASS <- trim.out$RD.PASS[inds]
  Seurat@meta.data$demux.RD.UNIQ <- trim.out$RD.UNIQ[inds]
  Seurat@meta.data$demux.N.SNP <- trim.out$N.SNP[inds]
  Seurat@meta.data$demux.PRB.DBL <- trim.out$PRB.DBL[inds]
  if(verbose){
    print("  Checking for barcode duplicates",quote=FALSE)
    print("",quote=FALSE)
    print("",quote=FALSE)}
  # Check if there are barcode duplicates accross 10X lanes.
  demux.barcode.dup <- array(FALSE, dim = length(cell.names))
  demux.barcode.dup[duplicated(inds, fromLast=FALSE)|duplicated(inds, fromLast=TRUE)] <- TRUE
  if (sum(demux.barcode.dup)>0){
    print("Warning: Cell barcodes were duplicated accross lanes of this dataset and may have been therefore artificially called as doublets by demuxlet. Recommended: Modify how your demuxlet was run to account for this.",quote=FALSE)
    Seurat@meta.data$demux.barcode.dup <- demux.barcode.dup
    print("",quote=FALSE)
    print("",quote=FALSE)
  }
  print("SUMMARY:",quote=FALSE)
  print(paste0(length(levels(lane.idents)),
               " lanes were identified and named:"),
        quote=FALSE)
  print(paste0(Lane.names[seq_along(levels(lane.idents))], collapse = ", "),
        quote=FALSE)
  print(paste0("The average number of SNPs per cell for all lanes was: ",round(mean(Seurat@meta.data$demux.N.SNP),1)),quote=FALSE)
  print(paste0("Out of ",length(cell.names)," cells total, Demuxlet assigned:"),quote=FALSE)
  print(paste0("     ",sum(Seurat@meta.data$demux.doublet.call=="SNG", na.rm = TRUE), " cells or ", round(100*sum(Seurat@meta.data$demux.doublet.call=="SNG", na.rm = TRUE)/length(cell.names),1), "% as singlets"),quote=FALSE)
  print(paste0("     ",sum(Seurat@meta.data$demux.doublet.call=="DBL", na.rm = TRUE), " cells or ", round(100*sum(Seurat@meta.data$demux.doublet.call=="DBL", na.rm = TRUE)/length(cell.names),1), "% as doublets"),quote=FALSE)
  print(paste0("     ",sum(Seurat@meta.data$demux.doublet.call=="AMB", na.rm = TRUE), " cells as too ambiguous to call."),quote=FALSE)
  print(paste0("     ",sum(!(cells %in% DEMUX.call$barcode)), " cells were not annotated in the demuxlet.best file."),quote=FALSE)
  Seurat
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
#' @export
demux.SNP.summary <- function(object = DEFAULT, group.by = "Lane", color.by = "Lane", plots = c("jitter","boxplot"), boxplot.color = "grey30", ...){
  DBPlot("demux.N.SNP", object, group.by = group.by, color.by = color.by, plots = plots, boxplot.color = boxplot.color, ...)
}

#' Plots the number of annotations per sample, per lane
#' @import ggplot2
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
#' @export
demux.calls.summary <- function(object = DEFAULT, singlets.only = TRUE,
                                main = "Sample Annotations by Lane", sub = NULL,
                                ylab = "Annotations", xlab = "Sample",
                                color = MYcolors[2], theme = NULL, rotate.labels = TRUE
){
  #Turn the object into a "name" if a full object was given
  if (typeof(object)=="S4") { object <- deparse(substitute(object)) }
  #Populate cells.use with a list of names, based on the singlets.only variable.
  cells.use <- if(singlets.only){
    .which_cells(meta("demux.doublet.call",object)=="SNG", object)
  } else {
    .all_cells(object)
  }
  #Establish the full list of cell/sample names
  all.cells <- .all_cells(object)
  #Set theme
  if(is.null(theme)){
    theme <- theme_bw()+theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              axis.line = element_blank(),
                              panel.border = element_blank())
  }

  dat <- as.data.frame.matrix(table(meta("Sample", object)[cells.use %in% all.cells],
                                    meta("Lane", object)[cells.use %in% all.cells]))
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
