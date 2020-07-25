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
            eval(expr(paste0(".data$",split.by[1], "~ .data$",split.by[2]))))
        )
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

.add_contours <- function(
    p, data, color, linetype = 1) {
    # Add contours based on the density of cells/samples
    # (Dim and Scatter plots)
    
    p + geom_density_2d(
        data = data,
        mapping = aes_string(x = "X", y = "Y"),
        color = color,
        linetype = linetype,
        na.rm = TRUE)
}

.add_labels <- function(
    p, Target_data, col.use = "color",
    labels.highlight, labels.size, labels.repel, split.by) {
    # Add text labels at/near the median x and y values for each group
    # (Dim and Scatter plots)

    #Determine medians
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
    }

    #Add labels
    args <- list(
        data = median.data,
        mapping = aes_string(x = "cent.x", y = "cent.y", label = "label"),
        size = labels.size)
    geom.use <-
        if (labels.highlight) {
            if (labels.repel) {
                ggrepel::geom_label_repel
            } else {
                geom_label
            }
        } else {
            if (labels.repel) {
                ggrepel::geom_text_repel
            } else {
                geom_text
            }
        }
    p + do.call(geom.use, args)
}

.calc_center_medians <- function(x.y.group.df, group.col) {
    groups <- levels(as.factor(as.character(x.y.group.df[,group.col])))
    data.frame(
        cent.x = vapply(
            groups,
            function(level) {
                median(x.y.group.df$X[x.y.group.df[,group.col]==level])
            }, FUN.VALUE = numeric(1)),
        cent.y = vapply(
            groups,
            function(level) {
                median(x.y.group.df$Y[x.y.group.df[,group.col]==level])
            }, FUN.VALUE = numeric(1)),
        label = groups)
}

.add_trajectory_lineages <- function(
    p, trajectories, clusters, arrow.size = 0.15, object, reduction.use,
    dim.1, dim.2) {
    # Add trajectory path arrows, following sets of cluster-to-cluster paths, from cluster median to cluster median.
    # (Dim and Scatter plots)
    #
    # p = a ggplot to add to
    # clusters = the name of the metadata slot that holds the clusters which were used for cluster-based trajectory analysis
    # trajectories = List of lists of cluster-to-cluster paths. Also, the output of SlingshotDataSet(SCE_with_slingshot)$lineages
    # arrow.size = numeric scalar that sets the arrow length (in inches) at the endpoints of trajectory lines.

    #Determine medians
    cluster.dat <- meta(clusters, object)
    cluster.levels <- metaLevels(clusters, object)
    data <- data.frame(
        cent.x = vapply(
            cluster.levels,
            function(level) {
                median(
                    .extract_Reduced_Dim(reduction.use, dim.1, object)$embedding[
                        cluster.dat==level])
            }, FUN.VALUE = numeric(1)),
        cent.y = vapply(
            cluster.levels,
            function(level) {
                median(
                    .extract_Reduced_Dim(reduction.use, dim.2, object)$embedding[
                        cluster.dat==level])
            }, FUN.VALUE = numeric(1)))

    #Add trajectories
    for (i in seq_along(trajectories)){
        p <- p + geom_path(
            data = data[as.character(trajectories[[i]]),],
            aes_string(x = "cent.x", y = "cent.y"),
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
                aes_string(x = "x", y = "y"),
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
                aes_string(x = "x", y = "y"),
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
            aes_string(x = "X", y = "Y", shape = col.use),
            color = "black", size=size*3/4, alpha = opacity) +
        scale_shape_manual(
            name = legend.title,
            values = letter.labels)
    p
}
