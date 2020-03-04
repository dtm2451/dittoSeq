#' Extracts Demuxlet information into a pre-made Seurat object
#' @importFrom utils read.table
#'
#' @param Seurat A pre-made Seurat object to add demuxlet information to.
#' @param lane.meta A meta.data
#' @param lane.names how the lanes should be named (if you want to give them something different from the default = Lane1, Lane2, Lane3...)
#' @param Demuxlet.best String pointing to the location of the .best output file from running of demuxlet OR a data.frame representing already imported .best matrix.
#'
#' Alternatively, can be a String vector representing the .best output locations of separate demuxlet runs which each correspond to 10X lanes that were combined together with cellranger aggr.  See the cellranger note below.
#' @param verbose whether to print messages about the stage of this process that is currently being run.
#' @param bypass.check for running this function anyway even when meta.data slots already around would be over-written.
#' @return A Seurat object with sample calls and relevant statistics added as metadata.
#' @details
#' The function takes in a previously generated Seurat object.
#' It also takes in demuxlet information either in the form of
#' 1: the location of a single demuxlet.best out file,
#' 2: the locations of multiple demuxlet.best output files,
#' or 3: a user-constructed data.frame created by reading in a demuxlet.best file.
#'
#' If a metadata slot name is provided to \code{lane.meta}, information in that metadata slot is copied into a metadata slot called "Lane".
#' Alternatively, if \code{lane.meta} is left as \code{NULL}, separate lanes are assumed to be marked by distinct values of "-#", as is the typical output of the 10X cellranger count & aggr pipeline.
#' In these situations, the \code{lane.names} input can be used to set specific names for each lane. "Lane1", "Lane2", "Lane3", etc, are used b y default.
#'
#' Barcodes (minus any "***_" often added by Seurat when objects are merged) are then matched between the \code{Seurat} and the \code{Demuxlet.best} dataframes.
#'
#' Note: "-#" information added by cellranger functions is not removed.
#' Doing so would cause cells, from separate 10X wells, which ended up with similar barcodes to become indistinguishable.
#' In demuxlet itself, ignorance of lane information leads to artificial doublet calls.
#' In \code{importDemux2Seurat}, ignorance of lane information can lead to import of improper demuxlet annotations.
#' For this reason, \code{importDemux2Seurat} checks for whether such artificial duplicates likely happened.
#' See the recommended cellranger/demuxlet pipeline below for specific suggestions for how to use this function with multi-well 10X data.
#'
#' Lane information and demuxlet calls and statistics are then imported into the \code{Seurat} as metadata:
#' \itemize{
#' \item Lane = guided by \code{lane.meta} import input, represents of separate droblet-generation lane, pool, sequencing lane, etc.
#' \item Sample = The sample call, parsed from the BEST column
#' \item demux.doublet.call = whether the sample was a singlet (SNG), doublet (DBL), or ambiguious (AMB), parsed from the BEST column
#' \item demux.RD.TOTL = RD.TOTL column
#' \item demux.RD.PASS = RD.PASS column
#' \item demux.RD.UNIQ = RD.UNIQ column
#' \item demux.N.SNP = N.SNP column
#' \item demux.PRB.DBL = PRB.DBL column
#' \item demux.barcode.dup = (Only generated when TRUEs will exist) whether a cell's barcode in the Demuxlet.best refered to only 1 cell in the Seurat object.
#' (When TRUE, indicates that cells from distinct lanes were interpretted together by demuxlet.
#' These will often be mistakenly called as doublets.)
#' }
#'
#' Finally, a summary of the results including mean number of SNPs and percentages of singlets and doublets is output unless \code{verbose} is set to \code{FALSE}.
#'
#' @note For use of multi-well 10X data + Demuxlet + this function:
#' 10X recommends running cellranger counts individually for each well/lane.
#' This leads to creation of separate genes x cells counts matrices for each lane.
#' *Demuxlet should also be run separately for each lane in order to minimize the informatic generation of artificial doublets.
#' Afterwards, there are two common methods of importing/merging such multi-well 10X data into Seurat.
#' Technical differences: Both options will alter the cell barcode names in a way that makes them unique across lanes, but how they do so is different.
#' Technical issue: Neither method adjusts the bacode names that are embedded within the BAM files which a user must supply to Demuxlet,
#' so that data needs to be modified in a proper way in order to make the Seurat cellnames and demuxlet BARCODEs match.
#'
#' \code{importDemux2Seurat} can work with either.
#' \itemize{
#' \item Option 1: merging matrices of all lanes with cellranger aggr before Seurat import.
#' Barcode uniquification method: A "-1", "-2", "-3", ... "-#" is appended to the end of all barcode names.
#' The number is incremented for each succesive lane. (Note: lane-numbers depend on the order in which they were supplied to cellranger aggr)
#' \item Option 2: Importing into Seurat, then merging Seurat objects.
#' Barcode uniquifiction method: user-defined strings are appended to the start of the barcodes, followed by an "_", within Seurat.
#' }
#'
#' Technical fix:
#' \code{importDemux2Seurat} ignores all information before a "_" in Seurat cellnames, but utilizes the "-#" information at the ends of Seurat cellnames.
#' \itemize{
#' \item Option 1: \code{importDemux2Seurat} can adjust the "-#" in the Demuxlet BARCODEs automatically for users before performing the matching step.
#' In order to take advantage of the automatic barcodes adjustment, just supply a vector containing the locations of the sepearate .best outputs for each lane, in the same order that lanes were combined in cellranger aggr.
#' \item Option 2: \code{importDemux2Seurat} ignores the user-adjustable lane information added by Seurat within cellnames.
#' To use with this method, run \code{importDemux2Seurat} on each lane's Seruat object separately & provide a unique name for each lane to the \code{lane.names} input, BEFORE merging into a single Seurat object.
#' }
#'
#' Run in these ways, demuxlet information can be matched to proper cells, and lane assignments can be properly reported in the "Lane" metadata slot.
#'
#' @seealso
#' Included QC visualizations:
#'
#' \code{\link{demux.calls.summary}} for plotting the number of sample annotations assigned within each lane.
#'
#' \code{\link{demux.SNP.summary}} for plotting the number of SNPs measured per cell.
#'
#' Or, see Kang et al. Nature Biotechnology, 2018. \url{https://www.nature.com/articles/nbt.4042}. For more information about the demuxlet cell-sample deconvolution method.
#' @examples
#'
#' pbmc <- Seurat::pbmc_small
#'
#' # Load in dittoSeq's included mock, (minimal) demuxlet output.
#' Demux <- dittoSeq::demuxlet.example
#'
#' ### Method 1 ###
#'
#' # Notice there is a groups metadata in this Seurat object.
#' getMetas(pbmc)
#' # We will treat these as if that holds Lane information
#'
#' # Now, running importDemux2Seurat:
#' pbmc <- importDemux2Seurat(
#'     Seurat = pbmc,
#'     lane.meta = "groups",
#'     Demuxlet.best = Demux)
#'
#' # Note, importDemux2Seurat can also take in the location of the .best file.
#' #   pbmc <- importDemux2Seurat(
#' #       Seurat = pbmc,
#' #       lane.meta = "groups",
#' #       Demuxlet.best = "Location/filename.best")
#'
#' # demux.SNP.summary() and demux.calls.summary() can now be used.
#' demux.SNP.summary(pbmc)
#' demux.calls.summary(pbmc)
#'
#' ### Method 2: cellranger aggr combined data ###
#'
#' # If cellranger aggr was used, lanes will be denoted by "-1", "-2", ... "-#"
#' #   at the ends of Seurat cellnames.
#' # Demuxlet should be run on each lane individually.
#' # Provided locations of each demuxlet.best output file, *in the same order
#' #   that lanes were provided to cellranger aggr* this function will then
#' #   adjust the "-#" within the .best BARCODEs automatically before matching
#' #
#' # pbmc <- importDemux2Seurat(
#' #     Seurat = pbmc,
#' #     Demuxlet.best = c(
#' #         "Location/filename1.best",
#' #         "Location/filename2.best"),
#' #     lane.names = c("g1","g2"))
#'
#' @author Daniel Bunis
#' @export

importDemux2Seurat <- function(
    Seurat, lane.meta = NULL, lane.names=NA, Demuxlet.best,
    bypass.check = FALSE, verbose = TRUE) {

    .error_if_no_Seurat()
    if (!verbose) {
        return(suppressMessages(importDemux2Seurat(
            Seurat,lane.meta,lane.names,Demuxlet.best,bypass.check)))
    }
    cell.names <- .all_cells(Seurat)

    .check_meta_overwrite(Seurat, bypass.check)

    message("Adding 'Lane' information as meta.data")
    Seurat <- .parse_and_add_lanes(Seurat, lane.meta, lane.names, cell.names)

    message("Extracting the Demuxlet calls")
    Demuxlet.info <- .extract_and_parse_demux_calls(Demuxlet.best)

    message("Matching barcodes")
    # Strip barcodes in cell.names from any of the extra info that may have been added by Seurat (normally "text_" at start of names)
    barcodes <- vapply(
        cell.names,
        function(X) strsplit(X, split = "_")[[1]][length(strsplit(X, split = "_")[[1]])],
        FUN.VALUE = character(1))
    # Remove any uneccesary cells from the Demux output (= cells not in this Seurat object)
    barcodes_in <- barcodes[barcodes %in% rownames(Demuxlet.info)]
    trim.info <- Demuxlet.info[barcodes_in,]
    # Determine the indices to match Seurat cells to Demux matrix.
    inds <- match(barcodes,barcodes_in)
    # Add an all NA column for any cells not found in the Demuxlet output
    trim.info <- rbind(trim.info, array(NA, dim = ncol(trim.info)))
    inds[is.na(inds)] <- nrow(trim.info)

    message("Adding Demuxlet info as metadata")
    Seurat <- .add_demux_metas_to_Seurat(Seurat, trim.info, inds)

    message("Checking for barcode duplicates across lanes...")
    Seurat <- .check_barcode_dups(Seurat, inds, NA.ind = nrow(trim.info))

    if (verbose) {
        .print_demux_import_summary(Seurat, cell.names, barcodes, barcodes_in)
    }

    Seurat
}

.check_meta_overwrite <- function(Seurat, bypass.check) {
    check.these <- c("Lane", "Sample", "demux.doublet.call", "demux.RD.TOTL",
        "demux.RD.PASS", "demux.RD.UNIQ", "demux.N.SNP", "demux.PRB.DBL",
        "demux.barcode.dup")
    if (sum(check.these %in% getMetas(Seurat))>0) {
        if(bypass.check) {
            message(
                "Note: '",
                paste0(
                    check.these[check.these %in% getMetas(Seurat)],
                    collapse = "' and '"),
                "' will be overwritten.\n")
        } else {
            stop(
                "WARNING: metadata slots would be overwritten:\n",
                paste0(
                    check.these[check.these %in% getMetas(Seurat)],
                    collapse = " and "),
                ".\n If this is okay, rerun with `bypass.check = TRUE`.")
        }
    }
}

.parse_and_add_lanes <- function(Seurat, lane.meta, lane.names, cell.names) {
    if (!(is.null(lane.meta))) {
        # Extract from given metadata
        lane.idents <- as.factor(
            as.character(meta(lane.meta, Seurat)))
        auto_lane.names <- metaLevels(lane.meta, Seurat)
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
    Seurat@meta.data$Lane <- lane.idents

    # Rename lanes from numbers to Lane# or given lane.names
    if (is.na(lane.names[1])) {
        lane.names <- auto_lane.names
    }
    levels(Seurat@meta.data$Lane) <- lane.names

    Seurat
}

.location_.best_to_data.frame <- function(locations) {
    read.demux <- function(file){
        message("    from \"", file, "\"")
        read.table(
            file = file,
            header=TRUE,
            sep="\t",
            stringsAsFactors = FALSE)
    }
    DF <- read.demux(locations[1])
    if (length(locations) > 1) {
        for (i in 2:length(locations)) {
            DF.new <- read.demux(locations[i])
            # Change cellranger count appended '-1' to cellranger aggr appended '-#'
            DF.new$BARCODE <- vapply(
                DF.new$BARCODE,
                function (barcode)
                    paste0(strsplit(barcode, split = "-1$")[[1]][1], "-", i),
                FUN.VALUE = character(1))
            DF <- cbind(DF, DF.new)
        }
    }

    DF
}

.extract_and_parse_demux_calls <- function(Demuxlet.best) {
    if (is.character(Demuxlet.best)) {
        Demuxlet.best <- .location_.best_to_data.frame(Demuxlet.best)
    }
    rownames(Demuxlet.best) <- Demuxlet.best$BARCODE
    Demuxlet.calls <- data.frame(t(
        vapply(
            as.character(Demuxlet.best$BEST),
            function(X) strsplit(X,'-')[[1]][seq_len(2)],
            FUN.VALUE = character(2))),
        row.names = as.character(Demuxlet.best$BARCODE),
        stringsAsFactors = FALSE)
    names(Demuxlet.calls) <- c("doublet", "sample")
    cbind(Demuxlet.best, Demuxlet.calls)
}

.add_demux_metas_to_Seurat <- function(Seurat, trim.info, inds) {
    Seurat@meta.data$Sample <- trim.info$sample[inds]
    Seurat@meta.data$demux.doublet.call <- trim.info$doublet[inds]
    Seurat@meta.data$demux.RD.TOTL <- trim.info$RD.TOTL[inds]
    Seurat@meta.data$demux.RD.PASS <- trim.info$RD.PASS[inds]
    Seurat@meta.data$demux.RD.UNIQ <- trim.info$RD.UNIQ[inds]
    Seurat@meta.data$demux.N.SNP <- trim.info$N.SNP[inds]
    Seurat@meta.data$demux.PRB.DBL <- trim.info$PRB.DBL[inds]
    Seurat
}

.check_barcode_dups <- function(Seurat, inds, NA.ind) {
    demux.barcode.dup <-
        duplicated(inds, fromLast=FALSE) | duplicated(inds, fromLast=TRUE)
    demux.barcode.dup[inds==NA.ind] <- NA
    if (sum(demux.barcode.dup, na.rm = TRUE)>0) {
        message("Warning: Cell barcodes were duplicated accross lanes of this",
            " dataset and may have been therefore artificially called as ",
            "doublets by demuxlet.")
        message("Recommended: Modify how demuxlet was run to account for this.\n")
        Seurat@meta.data$demux.barcode.dup <- demux.barcode.dup
    } else {
        message("  No barcode duplicates were found.\n")
    }
    Seurat
}

.print_demux_import_summary <- function(
    Seurat, cell.names, barcodes, barcodes_in) {

    num_singlets <- sum(meta("demux.doublet.call",Seurat)=="SNG", na.rm = TRUE)
    num_doublets <- sum(meta("demux.doublet.call",Seurat)=="DBL", na.rm = TRUE)
    num_ambig <- sum(meta("demux.doublet.call",Seurat)=="AMB", na.rm = TRUE)
    message("SUMMARY:\n",
        length(metaLevels("Lane", Seurat)),
        " lanes were identified and named:\n  ",

        paste0(metaLevels("Lane", Seurat), collapse = ", "),
        "\nThe average number of SNPs per cell for all lanes was: ",
        round(mean(meta("demux.N.SNP", Seurat), na.rm = TRUE),1), "\n",

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
#' Alternatively, a plotly object if \code{data.hover = TRUE} is provided, or list containing a ggplot and a dataframe if \code{data.out = TRUE} is provided.
#' @details
#' This function is a wrapper that essentially runs \code{\link{dittoPlot}}\code{("demux.N.SNP")} with a few modified defaults.
#' The defaults:
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
#' \code{\link{importDemux2Seurat}}, for how to import relevant demuxlet information as metadata.
#'
#' Kang et al. Nature Biotechnology, 2018. \url{https://www.nature.com/articles/nbt.4042}. For more information about the demuxlet cell-sample deconvolution method.
#' @examples
#' pbmc <- Seurat::pbmc_small
#'
#' # Generate some mock lane, sample, and doublet.call metadata
#' pbmc$demux.N.SNP <- rnorm(ncol(pbmc),75, 5)
#' pbmc$Lane <- sample(
#'     c("Lane1", "Lane2", "Lane3", "Lane4"),
#'     ncol(pbmc),
#'     replace = TRUE)
#'
#' demux.SNP.summary(pbmc)
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
#' @return For a given Seurat object, summarizes how many cells in each lane were anotated to each sample.  Assumes that the Sample calls of each cells, and which lane each cell belonged to, are stored in 'Sample' and 'Lane' metadata slots, respectively, as would be the case if the Seurat object was created with the Import10XDemux function.
#' @seealso
#' \code{\link{demux.SNP.summary}} for plotting the number of SNPs measured per cell.
#' This is the other Demuxlet-associated QC visualization included with dittoSeq.
#'
#' \code{\link{importDemux2Seurat}}, for how to import relevant demuxlet information as metadata.
#'
#' Kang et al. Nature Biotechnology, 2018. \url{https://www.nature.com/articles/nbt.4042}. For more information about the demuxlet cell-sample deconvolution method.
#' @examples
#' pbmc <- Seurat::pbmc_small
#'
#' # Generate some mock lane, sample, and doublet.call metadata
#' pbmc$Lane <- sample(
#'     c("Lane1", "Lane2", "Lane3", "Lane4"),
#'     ncol(pbmc),
#'     replace = TRUE)
#' pbmc$Sample <- sample(
#'     c("Sample-1", "Sample-2", "Sample-3", "Sample-4", "Sample-5"),
#'     ncol(pbmc),
#'     replace = TRUE)
#' pbmc$demux.doublet.call <- "SNG"
#'
#' demux.calls.summary(pbmc)
#'
#' demux.calls.summary(pbmc, data.out = TRUE)
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

    # Grab the data
    dat <- as.data.frame.matrix(
        table(
            as.character(meta("Sample", object)[cells.use %in% all.cells]),
            meta("Lane", object)[cells.use %in% all.cells]))
    dat$Sample <- row.names(dat)
    dat.m <- reshape2::melt(dat, "Sample")
    colnames(dat.m) <- c("Sample", "Lane", "Counts")

    if (data.out) {
        return(dat.m)
    } else {
        p <- ggplot(data = dat.m) + xlab(xlab) + ylab(ylab) + theme +
            ggtitle(main, subtitle = sub) + geom_hline(yintercept=0) +
            geom_col(
                aes_string(x = "Sample", y = "Counts"), fill = color) +
            facet_wrap(
                facets = ~Lane, ncol = 1, scales = "fixed", strip.position = "left")
        return(p)
    }
}
