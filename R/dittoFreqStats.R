#' Calculate per-sample frequencies of clusters or cell annotations, and compare them across group.
#' @inheritParams dittoPlot
#' @param cell.by String name of a per-cell metadata (a column of \code{colData(object)}) containing the cluster or cell annotation identities to quantify and assess.
#' @param group.by String name of a per-cell metadata (a column of \code{colData(object)}) containing sample-group identities.
#' @param group.1,group.2 Strings naming the 2 groups within the \code{group.by} metadata which you aim to compare.
#' @param sample.by String name of a per-cell metadata (a column of \code{colData(object)}) containing which sample each cell belongs to.
#' Recommendations for cyclone data, (standardiazed because they are required elements of the file_metadata input!): \itemize{
#' \item 'file_name': holds which original fcs file each cell came from.
#' \item 'donor_id': holds which patient/mouse each cell came from.
#' \item somthing else: sometimes, your data might both break up samples' data acquisition accross multiple .fcs files (so 'file_name' would then be too specific) & contain multiple timepoints or conditions per sample (so 'donor_id' is not specific enough).
#' In such a case, the burden lies on the user to create a viable metadata (\code{<object>$<metadata-name> <- <properly-uniqued-values>}, and then use \code{sample.by = <metadata-name>})
#' \item 'donor_id' & data subsetting with `cells.use`: As an alternative to creating a new metadata to use for \code{sample.by}, subsetting to only cells from a specific timepoint or condition might serve a dual purpose of achieving your specific analysis goal && allowing 'donor_id' to properly scope to individual samples.
#' See \code{cells.use} input description for further details.
#' }
#' @param cell.targs (Optional) Single string or a string vector naming which cell groups of the \code{cell.by} metadata which should be targetted.
#' When not provided, the function will loop through all cell groups in the \code{cell.by} metadata.
#' @param pseudocount Number, zero by default. A value to add, within fold_change calculations only, to both \code{group.1} and \code{group.2} median frequencies in order to avoid division by zero errors.
#' When needed, we recommend something small relative to the lowest expected cell frequencies of the data, 0.000001 perhaps.
#' Although a relatively small value like this can lead to heavily inflated log fold change values in the extreme cases where \code{group.1} or \code{group.2} frequencies are 0 or near 0, a tiny pseudocount leaves all other fold change values only minimally affected.
#' @param p.adjust.method String, passed along to the \code{method} input of \code{\link[stats]{p.adjust}}, any valid option for that input will work. "fdr" by default.
#' @param data.out Logical. When set to \code{TRUE}, changes the output from the stats data.frame alone to a named list containing both the stats ("stats") and the underlying per-sample frequency calculations ("data").
#' @return a data.frame, or if \code{data.out} was set to \code{TRUE}, a named list containing 2 data.frames, 'stats' and the underlying 'data'.
#' @details The function starts by utilizing \code{\link{dittoFreqPlot}} for
#' \code{cell.by}-cell frequency calculation within \code{sample.by}-samples,
#' percent normalization,
#' marking which \code{sample.by}-samples belong to which \code{group.by}-groups,
#' and trimming to only: 1. requested \code{cell.targs}, 2. \code{group.by}-groups \code{group.1} and \code{group.2}, and 3. cells matching the \code{cells.use} requirements if any were given.
#' It then removes some unnecessary columns from the data.frame directly returned by \code{\link{dittoFreqPlot}}. (Set \code{data.out = TRUE} to see what this cleaned return looks like!)
#'
#' Afterwards, it loops through all \code{cell.targs}, building a row of the eventual stats return for each.
#' Of note, a \code{pseudocount} can be introduced in median fold change calculation to prevent errors from division by zero. The use of a \code{pseudocount} has no effect on p-values.
#' Lastly, \code{p.adjust.method} correction, FDR by default, is applied to the 'p' column and added as a 'padj' column before data is returned.
#' @section The stats data.frame return:
#' Each row holds statistics for an individual comparison.
#' The columns represent:
#' \itemize{
#' \item cell_group: this row's cluster or cell-annotation
#' \item comparison: this groups of \code{group.by} compared in this row, formatted \code{<group.1>_vs_<group.2>}.
#' (For compatibility with running the function multiple times, each targwtting distinct groups, and then concatenating all outputs together!)
#' \item median_g1: the median frequency for the given cell_group within samples from \code{group.1}
#' \item median_g2: the median frequency for the given cell_group within samples from \code{group.2}
#' \item median_fold_change: \code{(median_g1 + pseudocount) / (median_g2 + pseudocount)}. Although zero by default, a small \code{pseudocount} can be set used to prevent division by zero error cases while only nominally affecting the value of other cases.
#' \item median_log2_fold_change: \code{log2( median_fold_change )}
#' \item positive_fc_means_up_in: Value = \code{group.1}, just a minor note to help remember the directionality of these fold changes!
#' \item p: The p-value associated with comparison of cell_group percent frequencies of group.1 samples versus group.2 samples using a Mann Whitney U Test / wilcoxon rank sum test (\code{\link[stats]{wilcox.test}}).
#' \item padj: p-values corrected by the chosen \code{p.adjust.method}, FDR by default, built from running \code{p.adjust(stats$p, method = p.adjust.method)} per all hypotheses tested in this call to the \code{freq_stats} function.
#' }
#' @author Daniel Bunis
#' @export
#' @importFrom stats wilcox.test
#' @importFrom stats p.adjust
#' @importFrom stats median
dittoFreqStats <- function(
    object,
    sample.by,
    cell.by,
    group.by, group.1, group.2,
    cell.targs = NULL,
    cells.use = TRUE,
    pseudocount = 0,
    p.adjust.method = "fdr",
    data.out = FALSE
) {
    
    if (is.null(cell.targs)) {
        cell.targs <- dittoSeq::metaLevels(cell.by, object)
    }
    
    # Collect stats with dittoSeq
    data <- dittoFreqPlot(
        object,
        var = cell.by,
        vars.use = cell.targs,
        sample.by = sample.by,
        group.by = group.by,
        cells.use = cells.use & object[[group.by, drop=TRUE]] %in% c(group.1, group.2),
        data.out = TRUE
    )$data
    
    # Clean
    data$var.data <- NULL # Column only needed for making the plot
    data$grouping <- NULL # Column only needed for making the plot, it's included in a column with the metadata's own name
    data$label.count.total.per.facet <- NULL # Not needed once used for percent calculation
    names(data)[1] <- "cell_group"
    
    # Here, we loop through all the cell_groups being targeted, 1- calculating stats and 2- building a data.frame during each iteration.
    #  The lapply call performs the iteration, and gathers the data.frames output by each iteration into a list.
    #  That list of data.frames created in our lapply is then 'rbind'ed into a single data.frame.
    stats <- do.call(
        rbind,
        lapply(
            unique(data$cell_group),
            function(clust) {
                data_use <- data[data$cell_group==clust,]
                g1s <- as.vector(data_use[[group.by]]==group.1)
                g2s <- as.vector(data_use[[group.by]]==group.2)
                new <- data.frame(
                    cell_group = clust,
                    comparison = paste0(group.1, "_vs_", group.2),
                    median_g1 = median(data_use$percent[g1s]),
                    median_g2 = median(data_use$percent[g2s]),
                    stringsAsFactors = FALSE
                )
                if (new$median_g2==0) {
                    warning("Looks like a 'pseudocount' will be needed to avoid division by zero errors. Try adding 'pseudocount = 0.000001' to your call and see the '?freq_stats' documentation for details.")
                }
                new$median_fold_change <- (new$median_g1 + pseudocount) / (new$median_g2 + pseudocount)
                new$median_log2_fold_change <- log2(new$median_fold_change)
                new$positive_fc_means_up_in <- group.1
                new$p <- wilcox.test(x=data_use$percent[g1s],
                                     y=data_use$percent[g2s])$p.value
                new
            })
    )
    
    # Apply FDR correction
    stats$padj <- p.adjust(stats$p, method = p.adjust.method)
    
    # Output
    if (data.out) {
        list(stats = stats, frequency_data = data)
    } else {
        stats
    }
}