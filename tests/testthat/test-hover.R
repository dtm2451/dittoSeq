# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-hover.R")

gene1 <- "gene1"
gene2 <- "gene2"
gene3 <- "gene3"
meta1 <- "score"
meta2 <- "clusters"

test_that("Showing hover.data works for ScatterPlot", {
    expect_s3_class(
        dittoScatterPlot(gene1, gene2, object = seurat, do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident")),
        "plotly")
})

test_that("Showing hover.data works for DimPlot", {
    expect_s3_class(
        dittoDimPlot(gene1, object = seurat, do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident")),
        "plotly")
    ### Manual Check: gene counts should become integers
    expect_s3_class(
        dittoDimPlot(gene1, object = seurat, do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident"),
            hover.slot = "counts"),
        "plotly")
})

test_that("Showing hover.data works for BarPlot", {
    expect_s3_class(
        dittoBarPlot(
            meta1, object = seurat,
            group.by = meta2,
            do.hover = TRUE),
        "plotly")
})

test_that("Showing hover.data works for Plot", {
    expect_s3_class(
        dittoPlot(
            gene1, object = seurat,
            group.by = meta2, color.by = meta2,
            do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident")),
        "plotly")
    expect_s3_class(
        dittoBoxPlot(
            gene1, object = seurat,
            group.by = meta2, color.by = meta2,
            do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident")),
        "plotly")
    expect_warning(
        dittoRidgePlot(
            gene1, object = seurat, plots = c("ridgeplot", "jitter"),
            group.by = meta2, color.by = meta2,
            do.hover = TRUE,
            hover.data = c(gene1,meta1,"ident")),
        NULL)
})

test_that("Showing hover.data works for VarsAcrossGroups", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(c(gene1, gene2, gene3),
            object = seurat, group.by = meta2,
            do.hover = TRUE),
        "plotly")
})
