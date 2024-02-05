#' Extracts Demuxlet information into a pre-made SingleCellExperiment or Seurat object
#' @importFrom utils read.table
#'
#' @param object A pre-made Seurat(v3+) or SingleCellExperiment object to add demuxlet information to.
#' @param raw.cell.names A string vector consisting of the raw cell barcodes of the object as they would have been output by cellranger aggr.
#' Format per cell.name = NNN...NNN-# where NNN...NNN are the cell barcode nucleotides, and # is the lane number.
#' This input should be used when additional information has been added directly into the cell names outside of Seurat's standard merge prefix: "user-text_".
#' @param lane.meta A string which names a metadata slot that contains which cells came from which droplet-generation wells.
#' @param lane.names String vector which sets how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param demuxlet.best String or String vector pointing to the location(s) of the .best output file from running of demuxlet.
#'
#' Alternatively, a data.frame representing an already imported .best matrix.
#' @param trim.before_ Logical which sets whether any characters in front of an "_" should be deleted from the \code{raw.cell.names} before matching with demuxlet barcodes.
#' @param verbose whether to print messages about the stage of this process that is currently being run & also the summary at the end.
#' @param bypass.check Logical which sets whether the function should run even when meta.data slots would be over-written.
#' @details
#' The function takes in a previously generated Seurat or SingleCellExperiment (SCE) object.
#' 
#' It also takes in demuxlet information either in the form of: 
#' (1) the location of a single demuxlet.best out file,
#' (2) the locations of multiple demuxlet.best output files,
#' (3) a user-constructed data.frame created by reading in a demuxlet.best file.
#'
#' Then it matches barcodes and adds demuxlet-information to the Seurat or SCE as metadata.
#' 
#' For a note on how best to utilize this function with multi-lane droplet-based data, see the devoted section below.
#' 
#' Specifically:
#' 
#' 1. If a metadata slot name is provided to \code{lane.meta}, information in that metadata slot is copied into a metadata slot called "Lane".
#' Alternatively, if \code{lane.meta} is left as \code{NULL}, separate lanes are assumed to be marked by distinct values of "-#" at the end of cell names, as is the typical output of the 10X cellranger count & aggr pipeline.
#' 
#' (1a. If \code{demuxlet.best} was provided as a set of separate file locations (recommended usage in conjunction with 'cellranger aggr'),
#' the "-#" at the ends of BARCODEs columns from these files are incremented on read-in so that they can match the incrementation applied by cellranger aggr.
#' See the section on multi-lane scRNAseq for more.)
#' 
#' 2. Barcodes in the demuxlet .best data are then matched to barcodes in the \code{object}.
#' The cell names, \code{colnames(object)}, are used by default for this matching,
#' but if these have been modified from what would have been given to demuxlet
#' -- outside of \code{-#} at the end or \code{***_}'s at the beginning, as can be added in common merge functions --
#' \code{raw.cell.names} can be provided and these cell names used instead.
#' 
#' 3. Singlet/doublet/ambiguous calls and sample identities (1st only for doublets) are parsed and carried into metadata.
#'
#' 4. Finally, a summary of the results including mean number of SNPs and percentages of singlets and doublets is output unless \code{verbose} is set to \code{FALSE}.
#'
#' @return The Seurat or SingleCellExperiment object with metadata added for "Sample" calls and other relevant statistics.
#' 
#' @section Metadata Added:
#' Lane information and demuxlet calls and statistics are imported into the \code{object} as these metadata:
#' \itemize{
#' \item Lane = guided by \code{lane.meta} import input or "-#"s in barcodes, represents the separate droplet-generation lanes.
#' \item Sample = The sample call, parsed from the BEST column
#' \item demux.doublet.call = whether the sample was a singlet (SNG), doublet (DBL), or ambiguious (AMB), parsed from the BEST column
#' \item demux.RD.TOTL = RD.TOTL column
#' \item demux.RD.PASS = RD.PASS column
#' \item demux.RD.UNIQ = RD.UNIQ column
#' \item demux.N.SNP = N.SNP column
#' \item demux.PRB.DBL = PRB.DBL column
#' \item demux.barcode.dup = (Only generated when TRUEs will exist) whether a cell's barcode in the demuxlet.best refered to only 1 cell in the \code{object}.
#' (When TRUE, indicates that cells from distinct lanes were interpretted together by demuxlet.
#' These will often be mistakenly called as doublets.)
#' }
#'
#' @section For data from multi-(droplet-gen-)lane scRNAseq:
#' There are many different ways such data might initially be processed which will affect its accessibility to \code{importDemux()}.
#' 
#' \strong{Initial Processing:}
#' 10X recommends running cellranger counts individually for each well/lane.
#' Non-10X droplet-based data from separate lanes should also be processed separately, at least for the steps of collecting reads for individual cells.
#' NOT processing such droplet lanes separately will create artificial doublets from cells that ended up with similar barcodes, but in separate droplet-gen lanes.
#' Thus, proper processing initially leads to creation of separate counts matrices for each droplet-generation lane.
#' 
#' \strong{Combining data from each lane:}
#' These per-lane counts matrices can be combined in various ways.
#' All options will alter the cell barcode names in a way that makes them unique across lanes, but this uniquification is achieved varies.
#' 
#' \emph{Counts table combination methods generally do not adjust adjust BAM files -- specifically the cell names embedded within the BAM files which is demuxlet uses for its BARCODEs column.
#' Thus cell names data may needs to be modified in a proper way in order to make the \code{object}'s cell names and \code{demuxlet.best}'s BARCODEs match.}
#' 
#' \strong{Running Demuxlet:}
#' Demuxlet should also be run, separately, on the BAM files of each individual lane.
#' Imporperly running demuxlet on a combined BAM file can lead to loss of lane information and then to generation of artificial doublet calls for cells of distinct wells that received simiar barcodes.
#' The BAM file associated with each demuxlet run is what is used for generating the BARCODE column of the demuxlet output.
#' 
#' \strong{How importDemux() handles barcode matching:}
#' \code{importDemux} is built to work with the 'cellranger aggr' pipeline by default, but can be used for demuxlet datasets processed differently as well (Option 2).
#' \itemize{
#' \item Option 1: When you merge matrices of all lanes with \strong{cellranger aggr before R import},
#' aggr's barcode uniquification method is to increment a "-1", "-2", "-3", ... "-#" that is appended to the end of all barcode names.
#' The number is incremented for each succesive lane. Note that lane-numbers depend on the order in which they were supplied to cellranger aggr.
#' \itemize{
#' \item \strong{to use:} Simply supply a \code{demuxlet.best} a vector containing the locations of the sepearate '.best' outputs for each lane, \emph{in the same order that lanes were provided to aggr}.
#' 
#' \code{importDemux} will adjust the "-#" in the \code{demuxlet.best} BARCODEs automatically before performing the matching step.
#' }
#' \item Option 2: When you instead \strong{import} your counts data into a Seurat or SingleCellExperiment, and \strong{then merge} the separate objects into one, the uniquifiction method is dependent on your particular method.
#' \itemize{
#' \item \strong{to use:} For these methods, it is easiest to
#' 1) \emph{import} your counts data,
#' 2) transfer in your demuxlet info with \emph{importDemux()} to each lane's object idividually (You can supply unique lane identifiers to the \code{lane.names} input.),
#' and then 3) \emph{merge} the separate objects.
#' }
#' \item Extra notes for any alternative cases:
#' \itemize{
#' \item For Seurat's \code{merge()}, user-defined strings can be appended to the start of the barcodes, followed by an "_".
#' By default, \code{importDemux()} will ignore these, but such ignorance can be controlled with the \code{trim.before_} input.
#' \item Alternatively, cell names that are consistent with the \code{demuxlet.best} BARCODEs can be supplied to the \code{raw.cell.names} input.
#' }
#' }
#'
#' @seealso
#' Included QC visualizations:
#'
#' \code{\link{demux.calls.summary}} for plotting the number of sample annotations assigned within each lane.
#'
#' \code{\link{demux.SNP.summary}} for plotting the number of SNPs measured per cell.
#'
#' Or, see Kang et al. Nature Biotechnology, 2018 \url{https://www.nature.com/articles/nbt.4042} for more information about the demuxlet cell-sample deconvolution method.
#' @examples
#'
#' #Prep: loading in an example dataset and sample demuxlet data
#' example("importDittoBulk", echo = FALSE)
#' demux <- demuxlet.example
#' colnames(myRNA) <- demux$BARCODE[seq_len(ncol(myRNA))]
#'
#' ###
#' ### Method 1: Lanes info stored in a metadata
#' ###
#'
#' # Notice there is a groups metadata in this Seurat object.
#' getMetas(myRNA)
#' # We will treat these as if that holds Lane information
#'
#' # Now, running importDemux:
#' myRNA <- importDemux(
#'     myRNA,
#'     lane.meta = "groups",
#'     demuxlet.best = demux)
#'
#' # Note, importDemux can also take in the location of the .best file.
#' #   myRNA <- importDemux(
#' #       object = myRNA,
#' #       lane.meta = "groups",
#' #       demuxlet.best = "Location/filename.best")
#'
#' # demux.SNP.summary() and demux.calls.summary() can now be used.
#' demux.SNP.summary(myRNA)
#' demux.calls.summary(myRNA)
#'
#' ###
#' ### Method 2: cellranger aggr combined data (denoted with "-#" in barcodes)
#' ###
#'
#' # If cellranger aggr was used, lanes will be denoted by "-1", "-2", ... "-#"
#' #   at the ends of Seurat cellnames.
#' # Demuxlet should be run on each lane individually.
#' # Provided locations of each demuxlet.best output file, *in the same order
#' #   that lanes were provided to cellranger aggr* this function will then
#' #   adjust the "-#" within the .best BARCODEs automatically before matching
#' #
#' # myRNA <- importDemux(
#' #     object = myRNA,
#' #     demuxlet.best = c(
#' #         "Location/filename1.best",
#' #         "Location/filename2.best"),
#' #     lane.names = c("g1","g2"))
#'
#' @author Daniel Bunis
#' @export

importDemux <- function(
    object, raw.cell.names = NULL, lane.meta = NULL, lane.names=NA,
    demuxlet.best, trim.before_ = TRUE, bypass.check = FALSE,
    verbose = TRUE) {

    if (is(object, "Seurat")) {
        .error_if_no_Seurat()
    }

    .check_meta_overwrite(object, bypass.check, verbose)

    .msg_if(verbose,"Adding 'Lane' information as meta.data")
    if (is.null(raw.cell.names)) {
        raw.cell.names <- .all_cells(object)
    }
    object <- .parse_and_add_lanes(object, lane.meta, lane.names, raw.cell.names)

    .msg_if(verbose,"Extracting the Demuxlet calls")
    Demuxlet.info <- .extract_and_parse_demux_calls(demuxlet.best, verbose)

    .msg_if(verbose,"Matching barcodes")

    # Strip barcodes in cell.names from any of the extra info that may have been added by Seurat (normally "text_" at start of names)
    barcodes <- raw.cell.names
    if (trim.before_) {
    barcodes <- vapply(
        raw.cell.names,
        function(X) strsplit(X, split = "_")[[1]][length(strsplit(X, split = "_")[[1]])],
        FUN.VALUE = character(1))
    }

    # Remove any cells from the Demux output that are not in this object
    barcodes_in <- barcodes[barcodes %in% rownames(Demuxlet.info)]
    if (length(barcodes_in)<1) {
        stop("No barcodes match between 'object' and 'demuxlet.best'")
    }
    trim.info <- Demuxlet.info[barcodes_in,]

    inds <- match(barcodes,barcodes_in)

    # Add an all NA column for any cells not found in the Demuxlet output
    trim.info <- rbind(trim.info, array(NA, dim = ncol(trim.info)))
    inds[is.na(inds)] <- nrow(trim.info)

    .msg_if(verbose,"Adding Demuxlet info as metadata")
    object <- .add_demux_metas(object, trim.info, inds)

    .msg_if(verbose,"Checking for barcode duplicates across lanes...")
    object <- .check_barcode_dups(object,inds, NA.ind=nrow(trim.info), verbose)

    if (verbose) {
        .print_demux_import_summary(
            object, raw.cell.names, barcodes, barcodes_in)
    }

    object
}

.check_meta_overwrite <- function(object, bypass.check, verbose) {
    check.these <- c("Lane", "Sample", "demux.doublet.call", "demux.RD.TOTL",
        "demux.RD.PASS", "demux.RD.UNIQ", "demux.N.SNP", "demux.PRB.DBL",
        "demux.barcode.dup")
    if (sum(check.these %in% getMetas(object))>0) {
        if(bypass.check) {
            .msg_if(verbose,
                "Note: '",
                paste0(
                    check.these[check.these %in% getMetas(object)],
                    collapse = "', '"),
                "' are being overwritten.\n")
        } else {
            stop(
                "metadata slots would be overwritten:\n",
                paste0(
                    check.these[check.these %in% getMetas(object)],
                    collapse = " and "),
                "\n If this is okay, rerun with 'bypass.check = TRUE'.")
        }
    }
}

.parse_and_add_lanes <- function(object, lane.meta, lane.names, cell.names) {

    # Parse cells' lane identities and default (auto) lane names
    if (!(is.null(lane.meta))) {
        # Extract from given metadata
        lane.idents <- as.factor(
            as.character(meta(lane.meta, object)))
        auto_lane.names <- metaLevels(lane.meta, object)
    } else {
        if (grepl("-",cell.names[1])) {
            # Use the `#` of `-#` parts of cellnames (cellranger aggr notation)
            lane.idents <- as.factor(
                vapply(
                    cell.names,
                    function(X){
                        chunks <- strsplit(X, split = "-")[[1]]
                        as.numeric(chunks[length(chunks)])
                    }, FUN.VALUE = numeric(1)))
        } else {
            # Treat as a single lane
            lane.idents <- as.factor(array(1, dim = length(cell.names)))
        }
        auto_lane.names <- paste0("Lane",seq_along(levels(lane.idents)))
    }

    # Add the Lane metadata, currently in number form
    object$Lane <- lane.idents

    # Rename lanes from numbers
    if (is.na(lane.names[1])) {
        lane.names <- auto_lane.names
    }
    levels(object$Lane) <- lane.names

    object
}

.location_.best_to_data.frame <- function(locations, verbose) {
    read.demux <- function(file, verbose){
        .msg_if(verbose,"    from \"", file, "\"")
        read.table(
            file = file,
            header=TRUE,
            sep="\t",
            stringsAsFactors = FALSE)
    }
    DF <- read.demux(locations[1], verbose)
    if (length(locations) > 1) {
        for (i in 2:length(locations)) {
            DF.new <- read.demux(locations[i], verbose)
            # Change cellranger count appended '-1' to cellranger aggr appended '-#'
            DF.new$BARCODE <- vapply(
                DF.new$BARCODE,
                function (barcode)
                    paste0(strsplit(barcode, split = "-1$")[[1]][1], "-", i),
                FUN.VALUE = character(1))
            DF <- rbind(DF, DF.new)
        }
    }

    DF
}

.extract_and_parse_demux_calls <- function(demuxlet.best, verbose) {
    if (is.character(demuxlet.best)) {
        demuxlet.best <- .location_.best_to_data.frame(demuxlet.best, verbose)
    }
    rownames(demuxlet.best) <- demuxlet.best$BARCODE
    Demuxlet.calls <- data.frame(t(
        vapply(
            as.character(demuxlet.best$BEST),
            function(X) strsplit(X,'-')[[1]][seq_len(2)],
            FUN.VALUE = character(2))),
        row.names = as.character(demuxlet.best$BARCODE),
        stringsAsFactors = FALSE)
    names(Demuxlet.calls) <- c("doublet", "sample")
    cbind(demuxlet.best, Demuxlet.calls)
}

.add_demux_metas <- function(object, trim.info, inds) {

    object$Sample <- trim.info$sample[inds]
    object$demux.doublet.call <- trim.info$doublet[inds]
    object$demux.RD.TOTL <- trim.info$RD.TOTL[inds]
    object$demux.RD.PASS <- trim.info$RD.PASS[inds]
    object$demux.RD.UNIQ <- trim.info$RD.UNIQ[inds]
    object$demux.N.SNP <- trim.info$N.SNP[inds]
    object$demux.PRB.DBL <- trim.info$PRB.DBL[inds]

    object
}

.check_barcode_dups <- function(object, inds, NA.ind, verbose) {
    demux.barcode.dup <-
        duplicated(inds, fromLast=FALSE) | duplicated(inds, fromLast=TRUE)
    demux.barcode.dup[inds==NA.ind] <- NA
    if (sum(demux.barcode.dup, na.rm = TRUE)>0) {
        warning("Warning: Cell barcodes are duplicated accross lanes, possibly leading to artificial doublet calls.")
        object$demux.barcode.dup <- demux.barcode.dup
    } else {
        .msg_if(verbose,"  No barcode duplicates were found.\n")
    }
    object
}

.print_demux_import_summary <- function(
    object, cell.names, barcodes, barcodes_in) {
    num_singlets <- sum(meta("demux.doublet.call",object)=="SNG", na.rm = TRUE)
    num_doublets <- sum(meta("demux.doublet.call",object)=="DBL", na.rm = TRUE)
    num_ambig <- sum(meta("demux.doublet.call",object)=="AMB", na.rm = TRUE)
    message("SUMMARY:\n",
        length(metaLevels("Lane", object)),
        " lanes were identified and named:\n  ",

        paste0(metaLevels("Lane", object), collapse = ", "),
        "\nThe average number of SNPs per cell for all lanes was: ",
        round(mean(meta("demux.N.SNP", object), na.rm = TRUE),1), "\n",

        "Out of ", length(cell.names),
        " cells in the Seurat object, Demuxlet assigned:\n",

        "    ", num_singlets, " cells or ",
        round(100*num_singlets/length(cell.names),1), "% as singlets\n",

        "    ", num_doublets, " cells or ",
        round(100*num_doublets/length(cell.names),1), "% as doublets\n",

        "    and ", num_ambig, " cells as too ambiguous to call.\n",

        sum(!(barcodes %in% barcodes_in)),
        " cells were not annotated in the demuxlet.best file.")
}

#' Plots the number of SNPs sequenced per droplet
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param group.by String "name" of a metadata to use for grouping values.
#' Default is "Lane".
#' @param color.by String "name" of a metadata to use for coloring.
#' Default is whatever was provided to \code{group.by}.
#' @param plots String vector which sets the types of plots to include: possibilities = "jitter", "boxplot", "vlnplot", "ridgeplot".
#' NOTE: The order matters, so use c("back","middle","front") when inputing multiple to put them in the order you want.
#' @param boxplot.color The color of the lines of the boxplot.
#' @param min numeric value which sets the minimum value shown on the y-axis.
#' @param add.line numeric value(s) where a dashed horizontal line should go.
#' Default = 50, a high confidence minimum number of SNPs per cell for highly accurate demuxlet sample deconvolution.
#' @param ... extra arguments passed to \code{\link{dittoPlot}}
#' @return A ggplot, made with \code{\link{dittoPlot}} showing a summary of how many SNPs were available to Demuxlet for each cell of a dataset.
#'
#' Alternatively, a plotly object if \code{data.hover = TRUE} is provided.
#'
#' Alternatively, list containing a ggplot and the underlying data as a dataframe if \code{data.out = TRUE} is provided.
#' @details
#' This function is a wrapper that essentially runs \code{\link{dittoPlot}}\code{("demux.N.SNP")} with a few modified defaults.
#' The altered defaults:
#' \itemize{
#' \item Data is grouped and colored by the "Lane" metadata (unless \code{group.by} or \code{color.by} are adjusted otherwise).
#' \item Data is displayed as boxplots with gray lines on top of dots for individual cells (unless \code{plots} or \code{boxplot.color} are adjusted otherwise).
#' \item The plot is set to have minimum y axis value of zero (unless \code{min} is adjusted otherwise).
#' \item A dashed line is added at the value 50, a very conservative minimum number of SNPs for high confidence sample calls (unless \code{add.line} is adjusted otherwise).
#' }
#' @seealso
#' \code{\link{demux.calls.summary}} for plotting the number of sample annotations assigned within each lane.
#' This is the other Demuxlet-associated QC visualization included with dittoSeq.
#'
#' \code{\link{dittoPlot}}, as \code{demux.SNP.summary} is essentially just a \code{dittoPlot} wrapper.
#'
#' \code{\link{importDemux}}, for how to import relevant demuxlet information as metadata.
#'
#' Kang et al. Nature Biotechnology, 2018 \url{https://www.nature.com/articles/nbt.4042} for more information about the demuxlet cell-sample deconvolution method.
#' @examples
#' example(importDemux, echo = FALSE)
#' demux.SNP.summary(myRNA)
#'
#' #Function wraps dittoPlot. See dittoPlot docs for more examples
#'
#' @author Daniel Bunis
#' @export
demux.SNP.summary <- function(
    object, group.by = "Lane", color.by = group.by,
    plots = c("jitter","boxplot"), boxplot.color = "grey30",
    add.line = 50, min = 0, ...) {

    dittoPlot(
        object, "demux.N.SNP", group.by, color.by,
        plots = plots, boxplot.color = boxplot.color,
        add.line = add.line, min = min, ...)
}

#' Plots the number of annotations per sample, per lane
#' @import ggplot2
#'
#' @param object A Seurat or SingleCellExperiment object
#' @param singlets.only Whether to only show data for cells called as singlets by demuxlet. Default is TRUE. Note: if doublets are included, only one of their sample calls will be used.
#' @param main plot title. Default = "Sample Annotations by Lane"
#' @param sub plot subtitle
#' @param ylab y axis label, default is "Annotations"
#' @param xlab x axis label, default is "Sample"
#' @param color bars color. Default is the dittoColors skyBlue.
#' @param theme A complete ggplot theme. Default is a slightly modified theme_bw().
#' @param rotate.labels whether sample names / x-axis labels should be rotated or not. Default is TRUE.
#' @param data.out Logical, whether underlying data for the plot should be output instead of the plot itself.
#' @return A faceted ggplot summarizing how many cells in each lane were anotated to each sample.
#' Assumes that the Sample calls of each cell, and which lane each cell belonged to, are stored in 'Sample' and 'Lane' metadata slots, respectively, as would be the case if demuxlet information was imported with \code{\link{importDemux}}.
#'
#' Alternatively, value will be a data.frame containing the underlying data if \code{data.out = TRUE} is provided.
#' @seealso
#' \code{\link{demux.SNP.summary}} for plotting the number of SNPs measured per cell.
#' This is the other Demuxlet-associated QC visualization included with dittoSeq.
#'
#' \code{\link{importDemux}}, for how to import relevant demuxlet information as metadata.
#'
#' Kang et al. Nature Biotechnology, 2018 \url{https://www.nature.com/articles/nbt.4042} for more information about the demuxlet cell-sample deconvolution method.
#' @examples
#' example(importDemux, echo = FALSE)
#'
#' demux.calls.summary(myRNA)
#'
#' # Exclude doublets by setting 'singlets only = TRUE'
#' demux.calls.summary(myRNA,
#'     singlets.only = TRUE)
#'
#' # To return the underlying data.frame
#' demux.calls.summary(myRNA, data.out = TRUE)
#'
#' @author Daniel Bunis
#' @export
demux.calls.summary <- function(
    object, singlets.only = FALSE,
    main = "Sample Annotations by Lane", sub = NULL, ylab = "Annotations",
    xlab = "Sample", color = dittoColors()[2], theme = NULL,
    rotate.labels = TRUE, data.out = FALSE) {

    if (singlets.only) {
        cells.use <- .which_cells(
            meta("demux.doublet.call", object)=="SNG", object)
    } else {
        cells.use <- .all_cells(object)
    }
    all.cells <- .all_cells(object)

    # Grab the data
    dat <- as.data.frame.matrix(
        table(
            as.character(meta("Sample", object)[all.cells %in% cells.use]),
            meta("Lane", object)[all.cells %in% cells.use]))
    dat$Sample <- row.names(dat)
    dat.m <- reshape2::melt(dat, "Sample")
    colnames(dat.m) <- c("Sample", "Lane", "Counts")

    if(is.null(theme)){
        x.args <- list(size=12)
        if (rotate.labels) {
            x.args <- c(x.args, list(angle=45, hjust = 1, vjust = 1))
        }
        theme <- theme_bw() +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_blank(),
                panel.border = element_blank(),
                axis.text.x= do.call(element_text, x.args))
    }

    if (data.out) {
        return(dat.m)
    } else {
        p <- ggplot(data = dat.m) + xlab(xlab) + ylab(ylab) + theme +
            ggtitle(main, subtitle = sub) + geom_hline(yintercept=0) +
            geom_col(
                aes(x = .data$Sample, y = .data$Counts), fill = color) +
            facet_wrap(
                facets = ~Lane, ncol = 1, scales = "fixed", strip.position = "left")
        return(p)
    }
}
