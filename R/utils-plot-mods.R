#' @importFrom stats as.formula
.add_splitting <- function(p, split.by, nrow, ncol, split.args,
                           enforce.grid=FALSE, grid.dir="col") {

    # Adds ggplot faceting to go with 'split.by' utilization.

    # When split.by is length 1, the shape is controlled with ncol & nrow
    if (length(split.by) == 1) {
        if (enforce.grid) {
            if (grid.dir=="col") {
                split.args$rows <- as.formula(paste(".~",split.by))
            } else {
                split.args$rows <- as.formula(paste(split.by,"~."))
            }
            return(p + do.call(facet_grid, split.args))
        }
        # When split.by is length 1, the shape is controlled with ncol & nrow
        split.args$facets <- split.by
        split.args$nrow <- nrow
        split.args$ncol <- ncol
        return(p + do.call(facet_wrap, split.args))
    }

    # When split.by is length 2, the first element is used for rows, and the
    # second element is used for columns.
    if (length(split.by) == 2) {
        split.args$rows <-
            as.formula(paste0(".data$", split.by[1], "~ .data$", split.by[2]))
        return(p + do.call(facet_grid, split.args))
    }
}

.remove_legend <- function(ggplot) {
    # Shorthand for ggplot legend removal
    ggplot + theme(legend.position = "none")
}

#' @importFrom cowplot ggdraw get_legend
.grab_legend <- function(ggplot) {
    # Obtains and plots just the legend of a ggplot
    cowplot::ggdraw(cowplot::get_legend(ggplot))
}
