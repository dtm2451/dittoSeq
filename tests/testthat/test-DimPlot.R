# Tests for dittoDimPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DimPlot.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
gene <- "gene1"
cont <- "number"
disc <- "groups"
disc2 <- "age"
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))
cols <- c("red", "blue", "yellow", "green", "black", "gray", "white")

test_that("dittoDimPlot can plot continuous or discrete data & raw or normalized expression", {
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            gene, object=seurat),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            gene, object=seurat,
            slot = "counts"),
        "ggplot")
})

test_that("dittoDimPlot (and reduction.use defaults) work for sce too", {
    # SCE
    expect_s3_class(
        dittoDimPlot(
            disc,  object=sce),
        "ggplot")
})

test_that("dittoDimPlot basic tweaks work", {
    # Manuel Check: big dots
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            size = 10),
        "ggplot")
    # Manuel Check: triangles
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            shape.panel = 17),
        "ggplot")
    # Manuel Check: see through large dots
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            size = 5,
            opacity = 0.5),
        "ggplot")
})

test_that("dittoDimPlot main legend can be removed or adjusted", {
    ### Manual Check: Legend removed
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            legend.show = FALSE),
        "ggplot")
    ### Manual Check: Legend title = "WOW"
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            legend.title = "WOW"),
        "ggplot")
    ### Manual Check: Legend symbols LARGE
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            legend.size = 15),
        "ggplot")
})

test_that("dittoDimPlots can be subset to show only certain cells/samples with any cells.use method", {
    expect_s3_class(
        c1 <- dittoDimPlot(
            disc, object=seurat,
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoDimPlot(
            disc, object=seurat,
            cells.use = cells.logical),
        "ggplot")
    expect_equal(c1,c2)
    c3 <- dittoDimPlot(
        disc, object=seurat,
        cells.use = 1:40,
        data.out = TRUE)
    expect_equal(c1,c2)
    expect_equal(c1,c3$p)
    expect_equal(nrow(c3$Target_data), 40)
    # And if we remove an entire grouping...
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            cells.use = meta(disc,seurat)!="A"),
        "ggplot")
})

test_that("dittoDimPlot shapes can be a metadata and the same as or distinct from var", {
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            shape.by = disc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            shape.by = disc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            shape.by = disc2),
        "ggplot")
})

test_that("dittoDimPlot shapes can be adjusted in many ways", {
    ### Manual check: Shapes should be triangle and diamond
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            shape.by = disc2, shape.panel= 17:19),
        "ggplot")
    ### Manual check: Shapes should be enlarged even more in the legend
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            shape.by = disc2, shape.legend.size = 10),
        "ggplot")
    ### Manual check: Shapes legend title should be removed
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            shape.by = disc2, shape.legend.title = NULL),
        "ggplot")
})

test_that("dittoDimPlot reduction.use can be changed", {
    ### Manuel Check: these should all look obviously distinct
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat, reduction.use = "PCA"),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat, reduction.use = "PCA",
            dim.1 = 3, dim.2 = 5),
        "ggplot")
})

test_that("dittoDimPlots colors can be adjusted for discrete data", {
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            color.panel = cols),
        "ggplot")
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            color.panel = cols[6:1]),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            color.panel = cols,
            colors = 6:1),
        "ggplot")
})

test_that("dittoDimPlots color scales can be adjusted for continuous data", {
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            min = -5, max = 100),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            legend.breaks = seq(10,60,10)),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            legend.breaks = seq(10,60,10),
            legend.breaks.labels = c("WOW",2:5,"HEY!")),
        "ggplot")
})

test_that("dittoDimPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: top and right plot outline removed
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            theme = theme_classic()),
        "ggplot")
})

test_that("dittoDimPlots discrete labels can be adjusted", {
    # Manual Check: 5:9
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            rename.var.groups = 5:9),
        "ggplot")
    # Manual Check: 3:6
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            shape.by = disc2, rename.shape.groups = 3:6),
        "ggplot")
})

test_that("dittoDimPlot can be labeled or circled", {
    ### Manual Check: Labels should repel in the first two (and move between
    # plots), and 1&3 with background, 2&4 without, 5: smaller labels
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.label = TRUE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.label = TRUE,
            labels.highlight = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.label = TRUE,
            labels.repel = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.label = TRUE,
            labels.highlight = FALSE,
            labels.repel = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.label = TRUE,
            labels.size = 3),
        "ggplot")
})

test_that("dittoDimPlot trajectory adding works", {
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            add.trajectory.lineages = list(
                c("B","A","C"),
                c("C","A")),
            trajectory.cluster.meta = disc,
            do.label = TRUE),
        "ggplot")
    # Manual Check: Arrows should move & GROW.
    expect_s3_class(
        dittoDimPlot(
            cont, object=seurat,
            add.trajectory.lineages = list(
                c("C","A")),
            trajectory.cluster.meta = disc,
            trajectory.arrow.size = 1),
        "ggplot")
    # Manual Check: Arrows should be detached from points
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10)
                )),
            trajectory.cluster.meta = disc),
        "ggplot")
})

test_that("dittoDimPlot lettering works", {
    ### Manual Check: Letters should be added
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.letter = TRUE, size = 3),
        "ggplot")
    ### Manual Check: see through dots and letters
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            do.letter = TRUE, size = 3,
            opacity = 0.5),
        "ggplot")
})

test_that("dittoDimPlot can remove axes numbers", {
    ### Manual Check: Numbers should be removed from the axes
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat, show.axes.numbers = FALSE),
        "ggplot")
})

test_that("dittoDimPlot plotting order can be ordered by the data", {
    out <- dittoDimPlot(disc, object=seurat, data.out = TRUE, size = 10, order = "decreasing")
    ### Manual Check: Orange aalways in front
    expect_s3_class(
        out$p,
        "ggplot")
    ### Manual Check: Dark blue aalways in front
    expect_equal(
        out$Target_data$color,
        rev(
            dittoDimPlot(
                disc, object=seurat, data.out = TRUE, size = 10, order = "increasing"
                )$Target_data$color)
    )
})

test_that("dittoDimPlot can add extra vars to dataframe", {
    df1 <- dittoDimPlot(
            disc, object=seurat,
            data.out = TRUE)[[2]]
    expect_s3_class(
        df2 <- dittoDimPlot(
            disc, object=seurat,
            extra.vars = c(gene, disc2), data.out = TRUE)[[2]],
        "data.frame")
    expect_equal(ncol(df1), 3)
    expect_equal(ncol(df2), 5)
})

test_that("dittoDimPlot genes can be different data types", {
    df <- dittoDimPlot(gene, object = seurat, data.out = TRUE,
        slot = "counts")
    expect_equal(
        df$Target_data$color,
        round(df$Target_data$color,0))
    df <- dittoDimPlot(gene, object = sce, data.out = TRUE,
        assay = "counts")
    expect_equal(
        df$Target_data$color,
        round(df$Target_data$color,0))
    df <- dittoDimPlot(gene, object = sce, data.out = TRUE,
        adjustment = "relative.to.max")
    expect_equal(
        0:1,
        range(df$Target_data$color))
})

test_that("dittoDimPlot adding contours", {
    expect_s3_class(dittoDimPlot(object=seurat, disc,
        do.contour = TRUE),
        "ggplot")
    
    ### Manual Check: Contour lines light blue and dashed
    expect_s3_class(dittoDimPlot(object=seurat, disc,
        do.contour = TRUE,
        contour.color = "lightblue", contour.linetype = "dashed"),
        "ggplot")
})

test_that("dittoDimPlot with and without rasterization produces identical plots", {
    ### Manual Check: Plots should appear identical
    expect_s3_class(dittoDimPlot(object=seurat, disc,
        do.raster = TRUE),
        "ggplot")
    
    expect_s3_class(dittoDimPlot(object=seurat, disc),
        "ggplot")
})

test_that("dittoDimPlot ignores do.letter/do.label/do.ellipse for continuous data", {
    expect_message(dittoDimPlot(object=seurat, cont,
        do.label = TRUE),
        "do.label was/were ignored for non-discrete data", fixed = TRUE)
    expect_message(dittoDimPlot(object=seurat, cont,
        do.letter = TRUE),
        "do.letter was/were ignored for non-discrete data", fixed = TRUE)
    expect_message(dittoDimPlot(object=seurat, cont,
        do.ellipse = TRUE),
        "do.ellipse was/were ignored for non-discrete data", fixed = TRUE)
    
    # No message for discrete data && MANUAAL CHECK: ellipse is drawn 
    expect_message(dittoDimPlot(object=seurat, disc,
        do.ellipse = TRUE),
        NA)
})

test_that("dittoDimPlot can be faceted with split.by (1 or 2 vars)", {
    # MANUAL CHECK: FACETING
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            split.by = disc2),
        "ggplot")
    # horizontal
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            split.nrow = 1),
        "ggplot")
    # vertical
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            split.ncol = 1),
        "ggplot")
    # Grid with rows=age, cols=groups
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc)),
        "ggplot")
    # Works with cells.use (should have grey cells)
    expect_s3_class(
        dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            cells.use = cells.logical),
        "ggplot")
})

test_that("dittoDimPlot added features work with single-metadata faceting", {
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            do.label = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            do.ellipse = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            do.letter = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            do.contour = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            add.trajectory.lineages = list(
                    c("C","A")),
            trajectory.cluster.meta = disc)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = disc2,
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10))))),
        NA)
})

test_that("dittoDimPlot added features work with double-metadata faceting", {
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            do.label = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            do.ellipse = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            do.letter = TRUE)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            add.trajectory.lineages = list(
                    c("C","A")),
            trajectory.cluster.meta = disc)),
        NA)
    expect_error(
        print(dittoDimPlot(
            disc, object=seurat,
            split.by = c(disc2,disc),
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10))))),
        NA)
})

test_that("dittoDimPlot swap.rownames works", {
    expect_s3_class(
        dittoDimPlot(sce, "gene1_symb", swap.rownames = "symbol"),
        "ggplot")
})

