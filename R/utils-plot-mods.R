.add_splitting <- function(p, split.by, nrow, ncol, object, cells.use) {
    if (length(split.by) == 1) {
        return(p + facet_wrap(split.by, nrow = nrow, ncol = ncol))
    }
    if (length(split.by) == 2) {
        return(p + facet_grid(
            rows = vars(meta(split.by[1], object)[cells.use]),
            cols = vars(meta(split.by[2], object)[cells.use])))
    }
}

.remove_legend <- function(ggplot) {
    ggplot + theme(legend.position = "none")
}

.grab_legend <- function(ggplot) {
    cowplot::ggdraw(cowplot::get_legend(ggplot))
}
