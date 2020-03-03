# Tests for dittoBarPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-BarPlot.R")

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
grp1 <- "orig.ident"
grp2 <- "RNA_snn_res.0.8"
grp3 <- "RNA_snn_res.1"
cells.names <- colnames(pbmc)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,40))
pbmc.bulk <- importDittoBulk(pbmc.se)

test_that("dittoBarPlot can quantify clustering of groupings in percent or raw count", {
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3, scale = "count"),
        "ggplot")
})

test_that("dittoBarPlot works for SCE", {
    # SCE
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc.se, group.by = grp3),
        "ggplot")
    # bulk SCE
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc.bulk, group.by = grp3),
        "ggplot")
})

test_that("dittoBarPlots can be subset to show only certain cells/samples with either cells.use method", {
    expect_s3_class(
        c1 <- dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            cells.use = cells.logical),
        "ggplot")
    expect_equal(c1,c2)
    # And if we remove an entire X grouping...
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            cells.use = meta(grp3,pbmc)!=0),
        "ggplot")
    # And if we remove an entire var grouping...
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            cells.use = meta(grp2,pbmc)!=0),
        "ggplot")
})

test_that("dittoBarPlot main legend can be removed", {
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            legend.show = FALSE),
        "ggplot")
})

test_that("dittoBarPlots colors can be adjusted", {
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            color.panel = dittoColors()[5:1]),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            colors = 5:1),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            color.panel = c("red","blue","yellow")),
        "ggplot")
})

test_that("dittoBarPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: plot should be boxed
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoBarPlots y-axis can be adjusted", {
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            y.breaks = seq(0,1,0.25)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            min = -0.5, max = 2,
            y.breaks = seq(0,1,0.25)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            scale = "count",
            y.breaks = seq(0,45,15)),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            scale = "count",
            min = -5, max = 55,
            y.breaks = seq(0,45,15)),
        "ggplot")
})

test_that("dittoBarPlot var-labels can be adjusted and reordered", {
    # Manual Check: groups changed to pikachu and libre.
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            var.labels.rename = c("pikachu", "libre")),
        "ggplot")
    # Manual Check: 1 on top of zero, with colors reversed too (orange still on top).
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            var.labels.reorder = 2:1),
        "ggplot")
})

test_that("dittoBarPlots x-labels can be adjusted", {
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            x.labels = 5:7),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            x.reorder = 3:1),
        "ggplot")
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: L->R, fewest(5), mid-count(6), high-count(7), with horizontal labels
    expect_s3_class(
        dittoBarPlot(
            grp2, pbmc, group.by = grp3,
            scale = "count",
            x.labels = 5:7, x.reorder = 3:1, x.labels.rotate = FALSE),
        "ggplot")
})
