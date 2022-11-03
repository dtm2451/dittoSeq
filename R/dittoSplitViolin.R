###/ Edits by Dan in roxygen chunk vs dittoPlot:
###/ removed dittoPlot params and replaced with an inherit
###/ filled in new params how I think new inputs might be able to work
###/ Added you as author plus note to remind us to cite SO authors / CC-SA-4.0 license.
#' Plots continuous data for customizeable cells'/samples' groupings on a y- (or x-) axis
#'
#' @inheritParams dittoPlot
#' @param var This input selects the primary data that will be displayed on the y-axis of the plot. Options: \Itemize{
#' \item When intending to compare between groups given with \code{split.vln.by}: Single string representing the name of a metadata or gene, OR a vector with length equal to the total number of cells/samples in the dataset.
#' Alternatively, a string vector naming multiple genes or metadata.
#' \item When intending to compare between 2 genes or metadata: Two strings (\code{c("target1","target2")}) naming the metadata/genes to use for the left and right sides of the split violin.
#' }
#' @param split.vln.by String (or NULL). The name of a metadata to use for setting groups of the split violin when comparing \code{var}-levels between groups.
#' Ignored when \code{length(var)==2}.
#' @param split.vln.1,split.vln.2 Single values of the metadata targeted by \code{split.vln.by} which sets the left (.1) and right (.2) groups to use in the split violin. If not given, the first two levels or unique values of the \code{split.vln.by}-data are used by default.
#' @param multivar.aes "split" or "group" (but unlike in other functions, cannot be "color"), the plot feature to utilize for displaying 'var' value when \code{var} is given multiple genes or metadata.
#' 
#' @return a ggplot where continuous data, grouped by sample, age, cluster, etc., shown on either the y-axis by a violin plot, boxplot, and/or jittered points, or on the x-axis by a ridgeplot with or without jittered points.
#'
#' Alternatively when \code{data.out=TRUE}, a list containing the plot ("p") and the underlying data as a dataframe ("data").
#'
#' Alternatively when \code{do.hover = TRUE}, a plotly converted version of the ggplot where additional data will be displayed when the cursor is hovered over jitter points.
#' @details
#' The function creates a dataframe containing the metadata or expression data associated with the given \code{var} (or if a vector of data is provided, that data).
#' On the discrete axis, data will be grouped by the metadata given to \code{group.by} and colored by the metadata given to \code{color.by}.
#' The \code{assay} and \code{slot} inputs can be used to change what expression data is used when displaying gene expression.
#' If a set of cells to use is indicated with the \code{cells.use} input, the data is subset to include only those cells before plotting.
#'
#' The \code{plots} argument determines the types of data representation that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' Inclusion of \code{"ridgeplot"} overrides \code{"boxplot"} and \code{"vlnplot"} presence and changes the plot to be horizontal.
#'
#' When \code{split.by} is provided the name of a metadata containing discrete data, separate plots will be produced representing each of the distinct groupings of the split.by data.
#'
#' \code{dittoRidgePlot}, \code{dittoRidgeJitter}, and \code{dittoBoxPlot} are included as wrappers of the basic \code{dittoPlot} function
#' that simply change the default for the \code{plots} input to be \code{"ridgeplot"}, \code{c("ridgeplot","jitter")}, or \code{c("boxplot","jitter")},
#' to make such plots even easier to produce.
#'
#' @section Many characteristics of the plot can be adjusted using discrete inputs:
#' The \code{plots} argument determines the types of \strong{data representation} that will be generated, as well as their order from back to front.
#' Options are \code{"jitter"}, \code{"boxplot"}, \code{"vlnplot"}, and \code{"ridgeplot"}.
#' 
#' Each plot type has specific associated options which are controlled by variables that start with their associated string.
#' For example, all jitter adjustments start with "\code{jitter.}", such as \code{jitter.size} and \code{jitter.width}.
#'
#' Inclusion of \code{"ridgeplot"} overrides \code{"boxplot"} and \code{"vlnplot"} presence and changes the plot to be horizontal.
#'
#' Additionally:
#'
#' \itemize{
#' \item \strong{Colors can be adjusted} with \code{color.panel}.
#' \item \strong{Subgroupings:} \code{color.by} can be utilized to split major \code{group.by} groupings into subgroups.
#' When this is done in y-axis plotting, dittoSeq automatically ensures the centers of all geoms will align,
#' but users will need to manually adjust \code{jitter.width} to less than 0.5/num_subgroups to avoid overlaps.
#' There are also three inputs through which one can use to control geom-center placement, but the easiest way to do all at once so is to just adjust \code{vlnplot.width}!
#' The other two: \code{boxplot.position.dodge}, and \code{jitter.position.dodge}.
#' \item \strong{Line(s) can be added} at single or multiple value(s) by providing these values to \code{add.line}.
#' Linetype and color are set with \code{line.linetype}, which is "dashed" by default, and \code{line.color}, which is "black" by default.
#' \item \strong{Titles and axes labels} can be adjusted with \code{main}, \code{sub}, \code{xlab}, \code{ylab}, and \code{legend.title} arguments.
#' \item The \strong{legend can be hidden} by setting \code{legend.show = FALSE}.
#' \item \strong{y-axis zoom and tick marks} can be adjusted using \code{min}, \code{max}, and \code{y.breaks}.
#' \item \strong{x-axis labels and groupings} can be changed / reordered using \code{x.labels} and \code{x.reorder}, and rotation of these labels can be turned on/off with \code{x.labels.rotate = TRUE/FALSE}.
#' \item \strong{Shapes used} in conjunction with \code{shape.by} can be adjusted with \code{shape.panel}.
#' \item Single or multiple \strong{additional per-cell features can be retrieved} and stashed within the underlying data using \code{extra.vars}.
#' This can be very useful for making manual additional alterations \emph{after} dittoSeq plot generation.
#' }
#' @seealso
#' \code{\link{multi_dittoPlot}} for easy creation of multiple dittoPlots each focusing on a different \code{var}.
#'
#' \code{\link{dittoPlotVarsAcrossGroups}} to create dittoPlots that show summarized expression (or values for metadata), accross groups, of multiple \code{vars} in a single plot.
#'
#' \code{\link{dittoRidgePlot}}, \code{\link{dittoRidgeJitter}}, and \code{\link{dittoBoxPlot}} for shortcuts to a few 'plots' input shortcuts
#'
#' @examples
#' example(importDittoBulk, echo = FALSE)
#' myRNA
#'
#' # Basic dittoplot, with jitter behind a vlnplot (looks better with more cells)
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint")
#'
#' # Color distinctly from the grouping variable using 'color.by'
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     color.by = "conditions")
#' dittoPlot(object = myRNA, var = "gene1", group.by = "conditions",
#'     color.by = "timepoint")
#'
#' # Update the 'plots' input to change / reorder the data representations
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"))
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("ridgeplot", "jitter"))
#'
#' ### Provided wrappers enable certain easy adjustments of the 'plots' parameter.
#' # Quickly make a Boxplot
#' dittoBoxPlot(myRNA, "gene1", group.by = "timepoint")
#' # Quickly make a Ridgeplot, with or without jitter
#' dittoRidgePlot(myRNA, "gene1", group.by = "timepoint")
#' dittoRidgeJitter(myRNA, "gene1", group.by = "timepoint")
#'
#' ### Additional Functionality
#' # Modify the look with intuitive inputs
#' dittoPlot(myRNA, "gene1", "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     boxplot.color = "white",
#'     main = "CD3E",
#'     legend.show = FALSE)
#'
#' # Data can also be split in other ways with 'shape.by' or 'split.by'
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     shape.by = "clustering",
#'     split.by = "SNP") # single split.by element
#' dittoPlot(object = myRNA, var = "gene1", group.by = "timepoint",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     split.by = c("groups","SNP")) # row and col split.by elements
#'
#' # Multiple genes or continuous metadata can also be plotted by giving them as
#' #   a vector to 'var'. One aesthetic of the plot will then be used to display
#' #   'var'-info, and you can control which (faceting / "split", x-axis grouping
#' #   / "group", or color / "color") with 'multivar.aes':
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"))
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"),
#'     multivar.aes = "group")
#' dittoPlot(object = myRNA, group.by = "timepoint",
#'     var = c("gene1", "gene2"),
#'     multivar.aes = "color")
#' 
#' # For faceting, instead of using 'split.by', the target data can alternatively
#' #   be given to 'extra.var' to have it added in the underlying dataframe, then
#' #   faceting can be added manually for extra flexibility
#' dittoPlot(myRNA, "gene1", "clustering",
#'     plots = c("vlnplot", "boxplot", "jitter"),
#'     extra.var = "SNP") + facet_wrap("SNP", ncol = 1, strip.position = "left")
#'
#' @authors Daniel Bunis and Rebecca Jaszczak (+ license-related text)
#' @export

###/ Edits by Dan in function vs dittoPlot:
###/ removed unneeded inputs and filled in new ones with proper defaults
###/ trimmed code, in top-level function only, that is certainly not needed
###/ removed "color" from multivar option as I'm 99% sure color will always be used for the split_violin & I'm not sure if the geom_split_violin is built to support a group aes (which is how a separate color.by could get implemented)?
###/ filled in new params' defaults as NULL because there are cases when they will not be needed.
###/ set new names for the internal functions
###/ Some additional thoughts that seem relevant:
###/   Perhaps split.vln.by does not need defaulting depending on final implementation.
###/     If not defaulted, we wouldn't need to write our own 'input is missing' error message for cases when it was needed but not given!
###/   UI conflict: To determine between using a var of length 2 as multivar versus to set the violin_split, can check either 'is.null(split.vln.by)' or 'exists(split.vln.by)'
dittoSplitViolin <- function(
    object,
    var,
    group.by,
    split.vln.by = NULL,
    split.vln.1 = NULL,
    split.violin.2 = NULL,
    split.by = NULL,
    extra.vars = NULL,
    cells.use = NULL,
    multivar.aes = c("split", "group"),
    multivar.split.dir = c("col", "row"),
    assay = .default_assay(object),
    slot = .default_slot(object),
    adjustment = NULL,
    swap.rownames = NULL,
    do.hover = FALSE,
    hover.data = var,
    color.panel = dittoColors(),
    colors = seq_along(color.panel),
    theme = theme_classic(),
    main = "make",
    sub = NULL,
    ylab = "make",
    y.breaks = NULL,
    min = NA,
    max = NA,
    xlab = "make",
    x.labels = NULL,
    x.labels.rotate = NA,
    x.reorder = NULL,
    split.nrow = NULL,
    split.ncol = NULL,
    split.adjust = list(),
    vlnplot.lineweight = 1,
    vlnplot.width = 1,
    vlnplot.scaling = "area",
    add.line = NULL,
    line.linetype = "dashed",
    line.color = "black",
    legend.show = TRUE,
    legend.title = "make",
    data.out = FALSE) {

    multivar.aes <- match.arg(multivar.aes)
    multivar.split.dir <- match.arg(multivar.split.dir)
    
    #Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    #Establish the full list of cell/sample names
    all.cells <- .all_cells(object)

    #Parse Title Defaults
    exp <- NULL
    if (isGene(var[1], object, assay)) {
        exp <- " expression"
    }
    ylab <- .leave_default_or_null(ylab,
        default = paste0(var,exp),
        null.if = !(length(var)==1 && is.character(var)))
    xlab <- .leave_default_or_null(xlab,
        default = group.by,
        null.if = multivar.aes=="group" && length(var)>1 && length(var) != length(all.cells))
    main <- .leave_default_or_null(main, var,
        null.if = !(length(var)==1 && is.character(var)))
    legend.title <- .leave_default_or_null(legend.title, var,
        null.if = is.null(shape.by))

    # Grab the data
    ###/ Edits needed: /###
    ###/ var can be two gene/metadata
    ###/ color.by not used
    ###/ split.vln.by needs pulling in
    ###/ split.vln.1 and .2 need to be filled in if not given, data subset to just these, and set them as factor levels
    ###/ do.hover and related inputs are not needed or at least are skippable/removable at this stage
    gather_out <- .dittoSplitViolin_data_gather(object, var, group.by, color.by,
        c(split.by,extra.vars), cells.use, assay, slot, adjustment,
        swap.rownames, do.hover, hover.data, x.reorder, x.labels,
        split.by, multivar.aes, multivar.split.dir)
    Target_data <- gather_out$Target_data
    split.by <- gather_out$split.by

    # Make the plot
    ###/ Probably even this chunk should go into the .ditto_split_violin plotter function below...
    ###/ it was set up this way for dittoPlot because these bits were universal whereas
    ###/ other bits were dependent on the intended plotting direction (vertical vs ridge/horizontal)
    p <- ggplot(Target_data, aes_string(fill="color")) +
        theme +
        scale_fill_manual(name = legend.title, values=color.panel[colors]) +
        ggtitle(main, sub)
    
    p <- .ditto_split_violin(
        p, Target_data, plots, xlab, ylab, shape.by, jitter.size,
        jitter.width, jitter.color, shape.panel, jitter.shape.legend.size,
        jitter.shape.legend.show, jitter.position.dodge,
        do.raster, raster.dpi,
        boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
        boxplot.position.dodge, boxplot.lineweight,
        vlnplot.lineweight, vlnplot.width, vlnplot.scaling,
        add.line, line.linetype, line.color,
        x.labels.rotate, do.hover, y.breaks, min, max, object)
    
    # Extra tweaks
    if (!is.null(split.by)) {
        p <- .add_splitting(
            p, split.by, split.nrow, split.ncol, split.adjust)
    }
    
    if (!legend.show) {
        p <- .remove_legend(p)
    }
    
    ###/ Probably won't be used? /###
    # if (do.hover) {
    #     p <- .warn_or_jitter_plotly(p, plots)
    # }
    
    # DONE. Return the plot +/- data
    if (data.out) {
        list(
            p = p,
            data = Target_data)
    } else {
        p
    }
}

.ditto_split_violin <- function(
    p, Target_data, plots, xlab, ylab, shape.by,
    jitter.size, jitter.width, jitter.color,shape.panel,
    jitter.shape.legend.size, jitter.shape.legend.show, jitter.position.dodge,
    do.raster, raster.dpi,
    boxplot.width, boxplot.color, boxplot.show.outliers, boxplot.fill,
    boxplot.position.dodge, boxplot.lineweight,
    vlnplot.lineweight, vlnplot.width, vlnplot.scaling, add.line,
    line.linetype, line.color, x.labels.rotate, do.hover, y.breaks, min, max,
    object) {
    # This function takes in a partial dittoPlot ggplot object without any data
    # overlay, and parses adding the main data visualizations.
    # Adds plots based on what is requested in plots, ordered by their order.

    # Now that we know the plot's direction, set direction & y-axis limits
    p <- p + aes_string(x = "grouping", y = "var.data")
    
    if (!is.null(y.breaks)) {
        p <- p + scale_y_continuous(breaks = y.breaks)
    }
    if (!is.na(min) || !is.na(max)) {
        p <- p + coord_cartesian(ylim=c(min,max))
    }

    # Add Plots
    for (i in seq_along(plots)) {
        if (plots[i] == "vlnplot") {
            p <- p + geom_violin(
                size = vlnplot.lineweight,
                width = vlnplot.width,
                scale = vlnplot.scaling,
                na.rm = TRUE)
        }

        if (plots[i] == "boxplot") {
            boxplot.args <- list(
                width = boxplot.width,
                color = boxplot.color,
                lwd = boxplot.lineweight,
                alpha = ifelse(boxplot.fill, 1, 0),
                position = position_dodge(width = boxplot.position.dodge),
                na.rm = TRUE)
            if (is.na(boxplot.show.outliers)) {
                boxplot.show.outliers <- ifelse("jitter" %in% plots, FALSE, TRUE)
            }
            if (!boxplot.show.outliers) {
                boxplot.args$outlier.shape <- NA
            }
            p <- p + do.call(geom_boxplot, boxplot.args)
        }

        if (plots[i] == "jitter") {
            
            # Create geom_jitter() arguments
            jitter.args <- list(
                position = position_jitterdodge(
                      jitter.width = jitter.width,
                      jitter.height = 0,
                      dodge.width = jitter.position.dodge,
                      seed = NA
                ),
                size=jitter.size,
                color = jitter.color)
            
            geom_for_jitter <- geom_jitter
            if (do.raster) {
                .error_if_no_ggrastr()
                geom_for_jitter <- ggrastr::geom_jitter_rast
                jitter.args$raster.dpi <- raster.dpi
            }
            
            jitter.aes.args <- list()
            if (do.hover) {
                jitter.aes.args$text <-  "hover.string"
            }
            
            #If shape.by metadata given, use it. Else, shapes[1] which = dots (16) by default
            if (!is.null(shape.by) && isMeta(shape.by, object)) {
                
                # Set shape in aes & set scales/theming.
                jitter.aes.args$shape <- shape.by
                
                p <- p + scale_shape_manual(
                    values = shape.panel[seq_along(metaLevels(shape.by, object, rownames(Target_data)))])
                
                if (!is.na(jitter.shape.legend.size)){
                    p <- p + guides(shape = guide_legend(
                        override.aes = list(size=jitter.shape.legend.size)))
                }
                if (jitter.shape.legend.show==FALSE){
                    p <- p + guides(shape = "none")
                }
                
            } else {
                # Set shape outside of aes
                jitter.args$shape <- shape.panel[1]
            }
            
            jitter.args$mapping <- do.call(aes_string, jitter.aes.args)
            
            if (do.hover) {
                p <- p + suppressWarnings(do.call(geom_for_jitter, jitter.args))
            } else {
                p <- p + do.call(geom_for_jitter, jitter.args)
            }
        }
    }

    # Add labels and, if requested, lines
    p <- p + xlab(xlab) + ylab(ylab)
    if (is.na(x.labels.rotate) || x.labels.rotate) {
        p <- p + theme(axis.text.x= element_text(angle=45, hjust = 1, vjust = 1))
    }
    if (!is.null(add.line)) {
        p <- p + geom_hline(yintercept=add.line, linetype= line.linetype, color = line.color)
    }

    p
}

###/ Edits needed: /###
###/ var can be two gene/metadata
###/ color.by not used (This aesthetic will ALWAYS be used for the split violin?)
###/ split.vln.by needs pulling in
###/ split.vln.1 and .2 need to be filled in if not given, data subset to just these, and set them as factor levels
###/ do.hover and related inputs are not needed or at least are skippable/removable at this stage
###/ Algorithm: \*Probably\* determine if  
.dittoSplitViolin_data_gather <- function(
    object, var, group.by, color.by,
    extra.vars, cells.use,
    assay, slot, adjustment,
    swap.rownames,
    do.hover, hover.data = c(var, extra.vars),
    x.reorder, x.labels,
    split.by, multivar.aes, multivar.split.dir) {

    all.cells <- .all_cells(object)
    
    # Support multiple genes/metadata
    # For multivar, not for splitting the violin by var
    if (length(var)>1 && length(var) != length(all.cells)) {
        Target_data <- do.call(rbind, lapply(
            var, function(this.var) {
                col <- switch(
                    multivar.aes, "split"="var", "group"="grouping", "color"="color")
                this.out <- .dittoSplitViolin_data_gather_inner(
                    object, this.var, group.by, color.by, extra.vars, cells.use,
                    all.cells, assay, slot, adjustment, swap.rownames, do.hover,
                    hover.data, x.reorder, x.labels
                )
                this.out[[col]] <- this.var
                this.out
            }
        ))
        if (multivar.aes == "split") {
            split.by <- .multivar_adjust_split_by(
                split.by, multivar.split.dir, multivar.col.name = "var")
        }
    } else {
        # Single var
        Target_data <- .dittoPlot_data_gather_inner(
            object, var, group.by, color.by, extra.vars, cells.use,
            all.cells, assay, slot, adjustment, swap.rownames, do.hover,
            hover.data, x.reorder, x.labels
        )
    }
    list(Target_data = Target_data, split.by = split.by)
}

.dittoSplitViolin_data_gather_inner <- function(
    object, var, group.by, color.by, extra.vars, cells.use, all.cells,
    assay, slot, adjustment, swap.rownames, do.hover, hover.data,
    x.reorder, x.labels
) {
    
    # Populate cells.use with a list of names if it was given anything else.
    cells.use <- .which_cells(cells.use, object)
    
    ### Make dataframe for storing the plotting data:
    full_data <- data.frame(
        var.data = .var_OR_get_meta_or_gene(
            var, object, assay, slot, adjustment, swap.rownames),
        grouping = meta(group.by, object),
        color = meta(color.by, object),
        row.names = all.cells)
    # Add split and extra data
    full_data <- .add_by_cell(full_data, extra.vars, extra.vars, object, assay,
        slot, adjustment, mult = TRUE)
    
    # Add hover strings
    if (do.hover) {
        full_data$hover.string <- .make_hover_strings_from_vars(
            hover.data, object, assay, slot, adjustment)
    }

    Target_data <- full_data[all.cells %in% cells.use,]
    # Reorder / Relabel grouping data
    Target_data$grouping <-
        .rename_and_or_reorder(Target_data$grouping, x.reorder, x.labels)
    
    Target_data
}

