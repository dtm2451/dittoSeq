# Tests for dittoScatterPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-ScatterPlot.R")

# Most ScatterPlot features are used/tested in the test-DimPlot, so this will look light.

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
gene <- "CD3E"
cont <- "number"
disc <- "RNA_snn_res.1"
disc2 <- "RNA_snn_res.0.8"
cells.names <- colnames(pbmc)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,40))
cols <- c("red", "blue", "yellow", "green", "black", "gray", "white")

test_that("dittoScatterPlot can plot genes or metadata", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, object = pbmc),
        "ggplot")
})

test_that("dittoScatterPlot can work with an sce", {
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, object = pbmc.se),
        "ggplot")
})

test_that("dittoScatterPlot can overlay colors, continuous or discrete", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, cont, object = pbmc),
        "ggplot")
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, disc, object = pbmc),
        "ggplot")
})

test_that("dittoScatterPlot can overlay shapes", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, NULL, disc, object = pbmc),
        "ggplot")
})

test_that("dittoScatterPlot gene display can utilize different data.types", {
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, gene, NULL, object = pbmc,
            data.type.x = "raw",
            data.type.y = "relative",
            data.type.color = "raw.normalized.to.max"),
        "ggplot")
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, gene, NULL, object = pbmc,
            do.hover = TRUE, hover.data = "PF4",
            data.type.x = "raw",
            data.type.y = "relative",
            data.type.color = "raw.normalized.to.max",
            hover.data.type = "scaled"),
        "plotly")
})

test_that("dittoScatterPlot gene display can utilize different data.types, with SingleCellExperiment objects", {
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, gene, NULL, object = pbmc.se,
            data.type.x = "raw",
            data.type.y = "relative",
            data.type.color = "normalized.to.max"),
        "ggplot")
})
