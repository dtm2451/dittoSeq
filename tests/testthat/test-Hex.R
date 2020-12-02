# Tests for dittoDimHex & dittoScatterHex functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Hex.R")

###########
# dittoDimHex relies on dittoScatterHex, so generally, we can test both via testing the Dim version
# When dittoDimHex works, Scatter can:
# - take in numeric vector data directly for 'x.var' & 'y.var'
# - output the plot and dataframe as list with 'data.out = TRUE'
###########

seurat$number <- as.numeric(seq_along(colnames(seurat)))
gene <- "gene1"
cont <- "number"
disc <- "groups"
disc2 <- "age"
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))
cols <- c("red", "blue", "yellow", "green", "black", "gray", "white")

test_that("DimHex & ScatterHex can plot density for Seurat or SCE", {
    expect_s3_class(dittoDimHex(object=seurat), "ggplot")
    expect_s3_class(dittoDimHex(object=sce), "ggplot")
    expect_s3_class(dittoDimHex(gene, cont, object=seurat), "ggplot")
    expect_s3_class(dittoDimHex(gene, cont, object=sce), "ggplot")
})

test_that("dittoDimHex - bins input adjusts number of bins", {
    ### Manual check: Large bins
    expect_s3_class(dittoDimHex(sce, bins = 5), "ggplot")
})

test_that("DimHex can plot continuous or discrete color.var data", {
    # Metadata
    expect_s3_class(dittoDimHex(object=seurat, disc), "ggplot")
    expect_s3_class(dittoDimHex(object=seurat, cont), "ggplot")
    # Expression
    expect_s3_class((p1 <- dittoDimHex(object=seurat, gene, data.out = TRUE))[[1]], "ggplot")
    expect_s3_class((p2 <- dittoDimHex(gene, object=seurat, slot = "counts", data.out = TRUE))$plot, "ggplot")
    expect_s3_class((p3 <- dittoDimHex(gene, object=sce, assay = "counts", data.out = TRUE))$plot, "ggplot")
    expect_s3_class((p4 <- dittoDimHex(gene, object=sce, adjustment = "relative.to.max", data.out = TRUE))$plot, "ggplot")
    expect_false(identical(p1$data$color, p2$data$color))
    expect_true(identical(p2$data$color, p3$data$color))
    expect_equal(max(p4$data$color), 1)
})

test_that("DimHex - color.method options work for discrete data, and defaults to 'max'", {
    ### Manual: Should have continuous color-scale and max.props in its title
    expect_s3_class(dittoDimHex(object=seurat, disc,
        color.method = "max.prop"),
        "ggplot")
    
    ### Manual: Next 2 should be the same plot with discrete color legend and "max" in its title
    expect_s3_class(dittoDimHex(object=seurat, disc,
        color.method = "max"),
        "ggplot")
    expect_s3_class(dittoDimHex(object=seurat, disc),
        "ggplot")
    
    expect_error(dittoDimHex(object=seurat, disc,
        color.method = "abcde"),
        "'color.method' not valid", fixed = TRUE)
})

test_that("DimHex - color.method options work for continuous data, and defaults to 'median'", {
    ### Manual: First should have lower max color than second
    expect_s3_class(dittoDimHex(object=seurat, cont,
        color.method = "max"),
        "ggplot")
    expect_s3_class(dittoDimHex(object=seurat, cont,
        color.method = "sum"),
        "ggplot")
    
    ### Manual: Next 2 should be the same plot
    expect_s3_class(dittoDimHex(object=seurat, cont,
        color.method = "median"),
        "ggplot")
    expect_s3_class(dittoDimHex(object=seurat, cont),
        "ggplot")
    
    expect_error(dittoDimHex(object=seurat, cont,
        color.method = "abcde"),
        "'color.method' not valid", fixed = TRUE)
})

test_that("DimHex dimensions can be adjusted", {
    expect_s3_class((p1 <- dittoDimHex(disc, object=seurat, data.out = TRUE,
        reduction.use = "PCA"))$plot,
        "ggplot")
    expect_s3_class((p2 <- dittoDimHex(disc, object=sce, data.out = TRUE,
        reduction.use = "PCA", dim.1 = 2, dim.2 = 1))$plot,
        "ggplot")
    expect_true(identical(p1$data$X, p2$data$Y))
    expect_true(identical(p1$data$Y, p2$data$X))
})

test_that("dittoDimHexs can be subset to show only certain cells/samples with any cells.use method", {
    expect_s3_class((c1 <- dittoDimHex(object=seurat, data.out = TRUE,
        cells.use = cells.names))$plot,
        "ggplot")
    expect_s3_class((c2 <- dittoDimHex(object=seurat, data.out = TRUE,
        cells.use = cells.logical))$plot,
        "ggplot")
    expect_s3_class((c3 <- dittoDimHex(object=seurat, data.out = TRUE,
        cells.use = 1:40))$plot,
        "ggplot")
    expect_equal(c1$data,c2$data)
    expect_equal(c1$data,c3$data)
    expect_equal(nrow(c3$data), 40)
    # And if we remove an entire grouping...
    expect_s3_class(dittoDimHex(disc, object=seurat,
            cells.use = meta(disc,seurat)!=0),
        "ggplot")
})

test_that("dittoDimHexs colors can be adjusted for discrete data", {
    expect_s3_class(dittoDimHex(object=seurat, disc,
        color.panel = cols), "ggplot")
    
    ### Manual check: These two should look the same.
    expect_s3_class(dittoDimHex(object=seurat, disc,
        color.panel = cols[5:1]), "ggplot")
    expect_s3_class(dittoDimHex(object=seurat, disc,
        color.panel = cols,
        colors = 5:1), "ggplot")
})

test_that("dittoDimHexs color legend: groupings can be renamed", {
    
    ### Manual check: color groups should be 1:5 (instead of A:E)
    expect_s3_class(dittoDimHex(object=seurat, disc,
        rename.color.groups = 1:5), "ggplot")
})

test_that("dittoDimHexs color scales can be adjusted for continuous color data", {
    
    ### Manual check: Legend range adjusted and black to orange
    expect_s3_class(dittoDimHex(object=seurat, cont,
        min = -5, max = 150, min.color = "black", max.color = "orange"),
        "ggplot")
    
    ### Manual check: Legend has breaaks at all tens in 10 to 60
    expect_s3_class(dittoDimHex(object=seurat, cont,
        legend.color.breaks = seq(10,60,10)),
        "ggplot")
    
    ### Manual check: Plot looks similar to above except from "WOW", 2:5, to "HEY"
    expect_s3_class(dittoDimHex(object=seurat, cont,
        legend.color.breaks = seq(10,60,10),
        legend.color.breaks.labels = c("WOW",2:5,"HEY!")),
        "ggplot")
})

test_that("dittoDimHexs color scales can be adjusted for density (color)", {
    
    ### Manual check: Legend range adjusted and black to orange
    expect_s3_class(dittoDimHex(object=seurat,
        min.density = -2, max.density = 2, min.color = "black", max.color = "orange"),
        "ggplot")
    
    ### Manual check: Legend from 1:3
    expect_s3_class(dittoDimHex(object=seurat,
        legend.density.breaks = seq(1:3)),
        "ggplot")
    
    ### Manual check: Plot looks similar to above except from "WOW", 2, to "HEY"
    expect_s3_class(dittoDimHex(object=seurat,
        legend.density.breaks = seq(1:3),
        legend.density.breaks.labels = c("WOW",2,"HEY!")),
        "ggplot")
})

test_that("dittoDimHexs color scales can be adjusted for density (opacity)", {
    
    ### Manual check: Opacity legend range adjusted -2 to 2 & barely any different
    expect_s3_class(dittoDimHex(object=seurat, cont,
        min.density = -2, max.density = 2, min.opacity = 0.5, max.opacity = 0.6),
        "ggplot")
    
    ### Manual check: Opacity legend breaks only at 1 and 3
    expect_s3_class(dittoDimHex(object=seurat, cont,
        legend.density.breaks = c(1,3)),
        "ggplot")
    
    ### Manual check: Opaacity legend from "WOW", 2, to "HEY"
    expect_s3_class(dittoDimHex(object=seurat, cont,
        legend.density.breaks = seq(1:3),
        legend.density.breaks.labels = c("WOW",2,"HEY!")),
        "ggplot")
})

test_that("dittoDimHexs titles and theme can be adjusted", {
    
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoDimHex(
            cont, object=seurat,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.color.title = "groups",
            legend.density.title = "Encounters"),
        "ggplot")
    
    ### Manual check: density legend (color)  = Encounters
    expect_s3_class(
        dittoDimHex(
            object=seurat,
            legend.density.title = "Encounters"),
        "ggplot")
    
    ### Manual check: top and right plot outline removed
    expect_s3_class(dittoDimHex(cont, object=seurat,
        theme = theme_classic()),
        "ggplot")
    
    ### Manual Check: Numbers should be removed from the axes
    expect_s3_class(dittoDimHex(disc, object=seurat,
        show.axes.numbers = FALSE),
        "ggplot")
    
    ### Manual Check: Legend removed
    expect_s3_class(dittoDimHex(object=seurat,
        legend.show = FALSE),
        "ggplot")
})

test_that("dittoDimHex can add extra vars to dataframe", {
    df1 <- dittoDimHex(object=seurat,
            data.out = TRUE)$data
    expect_s3_class(
        df2 <- dittoDimHex(object=seurat,
            extra.vars = c(gene, disc2), data.out = TRUE)$data,
        "data.frame")
    expect_equal(ncol(df1), 2)
    expect_equal(ncol(df2), 4)
})

test_that("dittoDimHex can be faceted with split.by (1 or 2 vars)", {
    
    # MANUAL CHECK: FACETING
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            split.by = disc2),
        "ggplot")
    
    # MANUAL CHECK: horizontal faceting
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            split.by = disc2,
            split.nrow = 1),
        "ggplot")
    
    # MANUAL CHECK: vertical faceting
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            split.by = disc2,
            split.ncol = 1),
        "ggplot")
    
    # MANUAL CHECK: Grid with rows=age, cols=groups
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            split.by = c(disc2,disc)),
        "ggplot")
    
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            split.by = c(disc2,disc),
            cells.use = cells.logical),
        "ggplot")
})

##########
# Added Features
##########

test_that("dittoDimHex trajectory adding works", {
    expect_s3_class(
        dittoDimHex(
            object=seurat, cont,
            add.trajectory.lineages = list(
                c("B","A","C"),
                c("C","A")),
            trajectory.cluster.meta = disc),
        "ggplot")
    
    ### Manual Check: One large arrow.
    expect_s3_class(
        dittoDimHex(
            object=seurat, cont,
            add.trajectory.lineages = list(
                c("C","A")),
            trajectory.cluster.meta = disc,
            trajectory.arrow.size = 1),
        "ggplot")
    
    ### Manual Check: Arrows should be detached from points
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10)
                ))),
        "ggplot")
})

test_that("dittoScatterHex trajectory curve adding works", {
    expect_s3_class(
        dittoScatterHex(
            gene, cont, gene, object=seurat,
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10)
                ))),
        "ggplot")
})

test_that("dittoDimHex adding contours", {
    expect_s3_class(dittoDimHex(object=seurat, disc,
        do.contour = TRUE),
        "ggplot")
    
    ### Manual Check: Contour lines light blue and dashed
    expect_s3_class(dittoDimHex(object=seurat, disc,
        do.contour = TRUE,
        contour.color = "lightblue", contour.linetype = "dashed"),
        "ggplot")
})

test_that("dittoDimHex do.label/do.ellipse", {
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            do.label = TRUE),
        "ggplot")
    expect_s3_class(
        dittoDimHex(
            disc, object=seurat,
            do.ellipse = TRUE),
        "ggplot")
})

test_that("dittoDimHex ignores do.label/do.ellipse for continuous data", {
    expect_message(dittoDimHex(object=seurat, cont,
        do.label = TRUE),
        "do.label was/were ignored for non-discrete data", fixed = TRUE)
    expect_message(dittoDimHex(object=seurat, cont,
        do.ellipse = TRUE),
        "do.ellipse was/were ignored for non-discrete data", fixed = TRUE)
    
    # No message for discrete data && MANUAAL CHECK: ellipse is drawn 
    expect_message(dittoDimHex(object=seurat, disc,
        do.ellipse = TRUE),
        NA)
})


##########
# Addition checks for Scatter
##########

# assay/slot/adjustment Scatter
test_that("dittoScatterHex gene display can utilize different data.types (excluding for hover)", {
    expect_s3_class((p <- dittoScatterHex(gene, gene, gene, object = seurat, data.out = TRUE,
        slot.x = "counts",
        slot.y = "counts",
        adjustment.color = "z-score"))$plot, "ggplot")
    expect_equal(
        p$data$X,
        round(p$data$Y,0))
    expect_equal(
        mean(p$data$color),
        0)
    expect_s3_class((p <- dittoScatterHex(gene, gene, gene, object = sce, data.out = TRUE,
        assay.x = "counts",
        assay.y = "counts"))$plot, "ggplot")
    expect_equal(
        p$data$X,
        round(p$data$Y,0))
    expect_s3_class((p <- dittoScatterHex(gene, gene, gene, object = sce, data.out = TRUE,
        adjustment.y= "relative.to.max"))$plot, "ggplot")
    expect_equal(
        max(p$data$Y), 1)
})

test_that("dittoScatterHex swap.rownames works", {
    expect_s3_class(
        dittoDimHex(sce, "gene1_symb", swap.rownames = "symbol"),
        "ggplot")
    expect_s3_class(
        dittoScatterHex(sce, "gene1_symb", "gene2_symb", "gene3_symb",
            swap.rownames = "symbol"),
        "ggplot")
})
