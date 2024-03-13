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

.add_letters_ellipses_labels_if_discrete <- function(
    p, data,
    is.discrete, do.letter, do.ellipse, do.label,
    labels.highlight, labels.size, labels.repel, labels.split.by,
    labels.repel.adjust,
    letter.size, letter.opacity, letter.legend.title, letter.legend.size,
    column = "color") {
    
    if (is.discrete) {

        if (do.letter) {
            p <- .add_letters(
                p, data, column,
                letter.size, letter.opacity, letter.legend.title, letter.legend.size)
        }
        
        if (do.ellipse) {
            p <- p + stat_ellipse(
                data=data,
                aes(x = .data$X, y = .data$Y, colour = .data[[column]]),
                type = "t", linetype = 2, linewidth = 0.5, show.legend = FALSE, na.rm = TRUE)
        }
        
        if (do.label) {
            p <- .add_labels(
                p, data, column, labels.highlight, labels.size,
                labels.repel, labels.repel.adjust, labels.split.by)
        }
        
    } else {
        
        # Data is incompatible, so message instead of adding.
        ignored.targs = paste(
            c("do.letter", "do.ellipse", "do.label")[c(do.letter,do.ellipse,do.label)],
            collapse = ", ")
        .msg_if(
            do.letter || do.ellipse || do.label,
            ignored.targs, " was/were ignored for non-discrete data.")
    }
    
    p
}

.add_contours <- function(
    p, data, color, linetype = 1) {
    # Add contours based on the density of cells/samples
    # (Dim and Scatter plots)
    
    p + geom_density_2d(
        data = data,
        mapping = aes(x = .data$X, y = .data$Y),
        color = color,
        linetype = linetype,
        na.rm = TRUE)
}

.add_labels <- function(
    p, Target_data, col.use = "color",
    labels.highlight, labels.size, labels.repel, labels.repel.adjust,
    split.by) {
    # Add text labels at/near the median x and y values for each group
    # (Dim and Scatter plots)

    # Determine medians
    if (is.null(split.by)) {
        
        median.data <- .calc_center_medians(Target_data, col.use)
        
    } else if (length(split.by)==1) {
        
        median.data <- NULL
        
        for (level in levels(as.factor(as.character(Target_data[,split.by])))) {
            
            level.dat <- Target_data[Target_data[,split.by]==level,]
            level.med.dat <- .calc_center_medians(level.dat, col.use)
            # Add split.by columns
            level.med.dat$split1 <- level
            colnames(level.med.dat)[4] <- split.by
            
            median.data <- rbind(median.data, level.med.dat)
        }
        
        # Ensure retention of factor level ordering
        median.data[,split.by] <- .retain_factor_level_order(
            median.data[,split.by], possible_factor = Target_data[,split.by])
        
    } else if (length(split.by)==2) {
        
        median.data <- NULL
        
        for (level1 in levels(as.factor(as.character(Target_data[,split.by[1]])))) {
            for (level2 in levels(as.factor(as.character(Target_data[,split.by[2]])))) {

                level.dat <- Target_data[Target_data[,split.by[1]]==level1,]
                level.dat <- level.dat[level.dat[,split.by[2]]==level2,]

                if (nrow(level.dat)>0) {
                    level.med.dat <- .calc_center_medians(level.dat, col.use)
                    # Add split.by columns
                    level.med.dat$split1 <- level1
                    level.med.dat$split2 <- level2
                    colnames(level.med.dat)[4:5] <- split.by
                    
                    median.data <- rbind(median.data, level.med.dat)
                }
            }
        }
        
        # Ensure retention of factor level ordering
        median.data[,split.by[1]] <- .retain_factor_level_order(
            median.data[,split.by[1]], possible_factor = Target_data[,split.by[1]])
        median.data[,split.by[2]] <- .retain_factor_level_order(
            median.data[,split.by[2]], possible_factor = Target_data[,split.by[2]])
    }

    #Add labels
    args <- list(
        data = median.data,
        mapping = aes(x = .data$cent.x, y = .data$cent.y, label = .data$label),
        size = labels.size)
    if (labels.repel) {
        if (is.list(labels.repel.adjust)) {
            args <- c(args, labels.repel.adjust)
        }
        geom.use <- if (labels.highlight) {
            ggrepel::geom_label_repel
        } else {
            ggrepel::geom_text_repel
        }
    } else {
        geom.use <- if (labels.highlight) {
            geom_label
        } else {
            geom_text
        }
    }

    p + do.call(geom.use, args)
}

.retain_factor_level_order <- function(new_data, possible_factor) {
    if (is.factor(possible_factor)) {
        factor(new_data, levels = levels(possible_factor))
    } else {
        new_data
    }
}

.calc_center_medians <- function(x.y.group.df, group.col) {
    groups <- levels(as.factor(as.character(x.y.group.df[,group.col])))
    data.frame(
        cent.x = vapply(
            groups,
            function(level) {
                median(x.y.group.df$X[x.y.group.df[,group.col]==level], na.rm = TRUE)
            }, FUN.VALUE = numeric(1)),
        cent.y = vapply(
            groups,
            function(level) {
                median(x.y.group.df$Y[x.y.group.df[,group.col]==level], na.rm = TRUE)
            }, FUN.VALUE = numeric(1)),
        label = groups)
}

.add_trajectory_lineages <- function(
    p, data, trajectories, clusters, arrow.size = 0.15, object) {
    # Add trajectory path arrows, following sets of cluster-to-cluster paths, from cluster median to cluster median.
    # (Dim and Scatter plots)
    #
    # p = a ggplot to add to
    # data = data for all cells/samples with columns including 'X' and 'Y'
    # clusters = the name of the metadata slot that holds the clusters which were used for cluster-based trajectory analysis.
    # trajectories = List of lists of cluster-to-cluster paths. Also, the output of SlingshotDataSet(SCE_with_slingshot)$lineages
    # arrow.size = numeric scalar that sets the arrow length (in inches) at the endpoints of trajectory lines.

    # Ensure data is in the objects' original cells' order, and add clusters data
    data <- data[.all_cells(object),]
    data$clusters <- meta(clusters, object)
    
    # Determine medians
    cluster.levels <- metaLevels(clusters, object)
    data <- .calc_center_medians(data, "clusters")

    #Add trajectories
    for (i in seq_along(trajectories)){
        p <- p + geom_path(
            data = data[as.character(trajectories[[i]]),],
            aes(x = .data$cent.x, y = .data$cent.y),
            arrow = arrow(
                angle = 20, type = "closed", length = unit(arrow.size, "inches")))
    }
    p
}

.add_trajectory_curves <- function(
    p, trajectories, arrow.size = 0.15, dim.1, dim.2) {
    # Add trajectory path arrows following sets of given (x,y) coordinates.
    # (Dim and Scatter plots)
    #
    # p = a ggplot to add to
    # trajectories = List of matrices containing trajectory curves. The output of SlingshotDataSet(SCE_with_slingshot)$curves can be used if the coordinate matrix (`$s`) for each list is extracted and they are all stored in a list.
    # arrow.size = numeric scalar that sets the arrow length (in inches) at the endpoints of trajectory lines.

    if ("s" %in% names(trajectories[[1]])) {
    #Add trajectories for princurves/slingshot list of lists provision method
        for (i in seq_along(trajectories)){
            #extract fit coords per cell
            data <- as.data.frame(trajectories[[i]]$s)
            #order cells' fit coords by pseudotime order
            data <- data[trajectories[[i]]$ord,]
            #name the dimensions used
            names(data)[c(dim.1,dim.2)] <- c("x", "y")
            p <- p + geom_path(
                data = data,
                aes(x = .data$x, y = .data$y),
                arrow = arrow(
                    angle = 20, type = "closed", length = unit(arrow.size, "inches")))
        }
    } else {
    #Add trajectories for general list of matrices provision method.
    #  Note: Accepts dataframes too.
        for (i in seq_along(trajectories)){
            data <- as.data.frame(trajectories[[i]])
            names(data) <- c("x", "y")
            p <- p + geom_path(
                data = data,
                aes(x = .data$x, y = .data$y),
                arrow = arrow(
                    angle = 20, type = "closed", length = unit(arrow.size, "inches")))
        }
    }
    p
}

.add_letters <- function(
    p, Target_data, col.use = "color", size, opacity, legend.title,
    legend.size) {
    # Overlay letters on top of the original colored dots.
    # Color blindness aid
    # (Dim and Scatter plots)

    letters.needed <- length(levels(as.factor(Target_data[,col.use])))
    letter.labels <- c(
        LETTERS, letters, 0:9, "!", "@", "#", "$", "%", "^", "&", "*", "(",
        ")", "-", "+", "_", "=", ";", "/", "|", "{", "}", "~"
        )[seq_len(letters.needed)]
    names(letter.labels) <- levels(as.factor(Target_data[,col.use]))
    p <- p +
        geom_point(
            data=Target_data,
            aes(x = .data$X, y = .data$Y, shape = .data[[col.use]]),
            color = "black", size=size*3/4, alpha = opacity) +
        scale_shape_manual(
            name = legend.title,
            values = letter.labels)
    p
}
