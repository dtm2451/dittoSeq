.add_splitting <- function(p, split.by, nrow, ncol, object, cells.use) {
    # Adds ggplot faceting to go with 'split.by' utilization.

    # When split.by is length 1, the shape is controlled with ncol & nrow
    if (length(split.by) == 1) {
        return(p + facet_wrap(split.by, nrow = nrow, ncol = ncol))
    }
    # When split.by is length 2, the first element is used for rows, and the
    # second element is used for columns.
    if (length(split.by) == 2) {
        return(p + facet_grid(
            rows = vars(meta(split.by[1], object)[cells.use]),
            cols = vars(meta(split.by[2], object)[cells.use])))
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
