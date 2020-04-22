# Tests for dittoBarPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-BarPlot.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
grp1 <- "groups"
grp2 <- "clusters"
grp3 <- "age"
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))

test_that("dittoBarPlot can quantify clustering of groupings in percent or raw count", {
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3, scale = "count"),
        "ggplot")
})

test_that("dittoBarPlot works for SCE", {
    # SCE
    expect_s3_class(
        dittoBarPlot(
            sce, grp2, group.by = grp3),
        "ggplot")
    # bulk SCE
    expect_s3_class(
        dittoBarPlot(
            bulk, grp2, group.by = grp3),
        "ggplot")
})

test_that("dittoBarPlots can be subset to show only certain cells/samples with any cells.use method", {
    expect_s3_class(
        c1 <- dittoBarPlot(
            seurat, grp2, group.by = grp3,
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoBarPlot(
            seurat, grp2, group.by = grp3,
            cells.use = cells.logical),
        "ggplot")
    expect_equal(c1,c2)
    expect_s3_class(
        c3 <- dittoBarPlot(
            seurat, grp2, group.by = grp3,
            cells.use = 1:40),
        "ggplot")
    expect_equal(c1,c3)
    # And if we remove an entire X grouping...
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            cells.use = meta(grp3,seurat)!=0),
        "ggplot")
    # And if we remove an entire var grouping...
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            cells.use = meta(grp2,seurat)!=0),
        "ggplot")
})

test_that("dittoBarPlot main legend can be removed", {
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            legend.show = FALSE),
        "ggplot")
})

test_that("dittoBarPlots colors can be adjusted", {
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            color.panel = dittoColors()[5:1]),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            colors = 5:1),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            color.panel = c("red","blue","yellow")),
        "ggplot")
})

test_that("dittoBarPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: plot should be boxed
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoBarPlots y-axis can be adjusted", {
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            y.breaks = seq(0,1,0.25)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            min = -0.5, max = 2,
            y.breaks = seq(0,1,0.25)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            scale = "count",
            y.breaks = seq(0,45,15)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            scale = "count",
            min = -5, max = 55,
            y.breaks = seq(0,45,15)),
        "ggplot")
})

test_that("dittoBarPlot var-labels can be adjusted and reordered", {
    # Manual Check: groups changed to pikachu and libre.
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            var.labels.rename = c("pikachu", "squirtle", "ivysaur", "charizard")),
        "ggplot")
    # Manual Check: 1 on top of zero, with colors reversed too (orange still on top).
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            var.labels.reorder = 4:1),
        "ggplot")
})

test_that("dittoBarPlots x-labels can be adjusted", {
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            x.labels = 4:7),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            x.reorder = 4:1),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: L->R 4:7, with horizontal labels
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            scale = "count",
            x.labels = 4:7, x.reorder = 4:1, x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: L->R 4:7, but order reversed compared to above
    expect_s3_class(
        dittoBarPlot(
            seurat, grp2, group.by = grp3,
            scale = "count",
            x.labels = 4:7, x.labels.rotate = FALSE),
        "ggplot")
})
