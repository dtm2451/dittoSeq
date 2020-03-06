# Tests for dittoScatterPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-ScatterPlot.R")

# Most ScatterPlot features are used/tested in the test-DimPlot, so this will look light.

seurat$number <- as.numeric(seq_along(colnames(seurat)))
gene <- "gene1"
cont <- "number"
disc <- "groups"
disc2 <- "age"
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))
cols <- c("red", "blue", "yellow", "green", "black", "gray", "white")

test_that("dittoScatterPlot can plot genes or metadata & work for SCE", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, object = seurat),
        "ggplot")
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, object = sce),
        "ggplot")
})

test_that("dittoScatterPlot can overlay colors, continuous or discrete", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, cont, object = seurat),
        "ggplot")
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, disc, object = seurat),
        "ggplot")
})

test_that("dittoScatterPlot can overlay shapes", {
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, NULL, disc, object = seurat),
        "ggplot")
})

test_that("dittoScatterPlot (and dittoDimPlot) can add extra vars to dataframe", {
    df1 <- dittoScatterPlot(
            gene, cont, NULL, disc, object = seurat,
            data.out = TRUE)[[2]]
    expect_s3_class(
        df2 <- dittoScatterPlot(
            gene, cont, NULL, disc, object = seurat,
            extra.vars = c(gene, disc2), data.out = TRUE)[[2]],
        "data.frame")
    expect_equal(ncol(df1), 3)
    expect_equal(ncol(df2), 5)
})

test_that("dittoScatterPlot (and dittoDimPlot) can be facted with split.by", {
    # MANUAL CHECK: FACETING
    expect_s3_class(
        dittoScatterPlot(
            gene, cont, NULL, disc, object = seurat,
            split.by = disc2),
        "ggplot")
})

test_that("dittoScatterPlot gene display can utilize different data.types (excluding for hover)", {
    expect_s3_class(
        dittoScatterPlot(
            gene, gene, gene, NULL, object = seurat,
            slot.x = "counts",
            adjustment.y = "z-score",
            adjustment.color = "relative.to.max"),
        "ggplot")
})
