# Tests for dittoDimPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DimPlot.R")

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
gene <- "CD3E"
cont <- "number"
disc <- "RNA_snn_res.1"
disc2 <- "RNA_snn_res.0.8"
cells.names <- colnames(pbmc)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,40))
cols <- c("red", "blue", "yellow", "green", "black", "gray", "white")

test_that("dittoDimPlot can plot continuous or discrete data & raw or normalized expression", {
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            gene, pbmc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            gene, pbmc,
            data.type = "raw"),
        "ggplot")
})

test_that("dittoDimPlot (and reduction.use defaults) work for sce too", {
    # SCE
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc.se),
        "ggplot")
    # RNAseq
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc.rnaseq),
        "ggplot")
})

test_that("dittoDimPlot basic tweaks work", {
    # Manuel Check: big dots
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            size = 10),
        "ggplot")
    # Manuel Check: triangles
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            shape.panel = 17),
        "ggplot")
    # Manuel Check: see through large dots
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            size = 5,
            opacity = 0.5),
        "ggplot")
})

test_that("dittoDimPlot main legend can be removed or adjusted", {
    ### Manual Check: Legend removed
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            legend.show = FALSE),
        "ggplot")
    ### Manual Check: Legend title = "WOW"
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            legend.title = "WOW"),
        "ggplot")
    ### Manual Check: Legend symbols LARGE
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            legend.size = 15),
        "ggplot")
})

test_that("dittoDimPlots can be subset to show only certain cells/samples with either cells.use method", {
    expect_s3_class(
        c1 <- dittoDimPlot(
            disc, pbmc,
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoDimPlot(
            disc, pbmc,
            cells.use = cells.logical),
        "ggplot")
    expect_equal(c1,c2)
    # And if we remove an entire grouping...
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            cells.use = meta(disc,pbmc)!=0),
        "ggplot")
})

test_that("dittoDimPlot shapes can be a metadata and the same as or distinct from var", {
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            shape.var = disc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            shape.var = disc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            shape.var = disc2),
        "ggplot")
})

test_that("dittoDimPlot shapes can be adjusted in many ways", {
    ### Manual check: Shapes should be triangle and diamond
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            shape.var = disc2, shape.panel= 17:19),
        "ggplot")
    ### Manual check: Shapes should be enlarged even more in the legend
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            shape.var = disc2, shape.legend.size = 10),
        "ggplot")
    ### Manual check: Shapes legend title should be removed
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            shape.var = disc2, shape.legend.title = NULL),
        "ggplot")
})

test_that("dittoDimPlot reduction.use can be changed", {
    ### Manuel Check: these should all look obviously distinct
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc, reduction.use = "pca"),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc, reduction.use = "pca",
            dim.1 = 3, dim.2 = 5),
        "ggplot")
})

test_that("dittoDimPlots colors can be adjusted for discrete data", {
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            color.panel = cols),
        "ggplot")
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            color.panel = cols[6:1]),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            color.panel = cols,
            colors = 6:1),
        "ggplot")
})

test_that("dittoDimPlots color scales can be adjusted for continuous data", {
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            min = -5, max = 100),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            legend.breaks = seq(10,60,10)),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            legend.breaks = seq(10,60,10),
            legend.breaks.labels = c("WOW",2:5,"HEY!")),
        "ggplot")
})

test_that("dittoDimPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: gridlines should be shown
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoDimPlots discrete labels can be adjusted", {
    # Manual Check: 5:7
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            rename.var.groups = 5:7),
        "ggplot")
    # Manual Check: 0:4
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            shape.var = disc2, rename.shape.groups = 3:4),
        "ggplot")
})

test_that("dittoDimPlot can be labeled or circled", {
    ### Manual Check: Labels should repel in the first two (and move between
    # plots), and 1&3 with background, 2&4 without, 5: smaller labels
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.label = TRUE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.label = TRUE,
            labels.highlight = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.label = TRUE,
            labels.repel = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.label = TRUE,
            labels.highlight = FALSE,
            labels.repel = FALSE),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.label = TRUE,
            labels.size = 3),
        "ggplot")
})

test_that("dittoDimPlot trajectory adding works", {
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            add.trajectory.lineages = list(
                c(1,0,2),
                c(2,1)),
            trajectory.cluster.meta = disc,
            do.label = TRUE),
        "ggplot")
    # Manual Check: Arrows should move & GROW.
    expect_s3_class(
        dittoDimPlot(
            cont, pbmc,
            add.trajectory.lineages = list(
                c(1,0)),
            trajectory.cluster.meta = disc,
            trajectory.arrow.size = 1),
        "ggplot")
    # Manual Check: Arrows should be detached from points
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            add.trajectory.curves = list(
                data.frame(
                    c(-10,0,-20),
                    c(-20,-10,0)),
                data.frame(
                    c(5:20),
                    c(5:10,9:5,6:10)
                )),
            trajectory.cluster.meta = disc,
            do.label = TRUE),
        "ggplot")
})

test_that("dittoDimPlot lettering works", {
    ### Manual Check: Letters should be added
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.letter = TRUE, size = 3),
        "ggplot")
    ### Manual Check: see through dots and letters
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc,
            do.letter = TRUE, size = 3,
            opacity = 0.5),
        "ggplot")
})

test_that("dittoDimPlot can remove axes numbers", {
    ### Manual Check: Numbers should be removed from the axes
    expect_s3_class(
        dittoDimPlot(
            disc, pbmc, show.axes.numbers = FALSE),
        "ggplot")
})
