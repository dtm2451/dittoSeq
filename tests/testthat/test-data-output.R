# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-data-output.R")

test_that("Data outputing works for ScatterPlot", {
    expect_type(
        dittoScatterPlot("gene1", "gene2", object = seurat,
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoScatterPlot("gene1", "gene2", object = seurat,
            data.out = TRUE)$Target_data,
        "list")
    expect_true(ncol(d1) >= 2 && nrow(d1) > 10)
})

test_that("Data outputing works for DimPlot", {
    expect_type(
        dittoDimPlot("gene1", object = seurat,
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoDimPlot("gene1", object = seurat,
            data.out = TRUE)$Target_data,
        "list")
    expect_true(ncol(d1) > 2 && nrow(d1) > 10)
})

test_that("Data outputing works for BarPlot", {
    expect_s3_class(
        dittoBarPlot(
            "clusters", object = seurat,
            group.by = "age", data.out = TRUE,
            do.hover = TRUE),
        "data.frame")
})

test_that("Data outputing works for Plot", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_type(
        dittoPlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoPlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE)$data,
        "list")
    expect_type(
        dittoRidgePlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE),
        "list")
    expect_type(
        d2 <- dittoRidgePlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE)$data,
        "list")
    expect_type(
        dittoBoxPlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE),
        "list")
    expect_type(
        d3 <- dittoBoxPlot(
            "gene1", object = seurat,
            group.by = "age", color.by = "age",
            data.out = TRUE)$data,
        "list")
    expect_true(ncol(d1) > 2 && nrow(d1) > 10)
    expect_true(ncol(d2) > 2 && nrow(d1) > 10)
    expect_true(ncol(d3) > 2 && nrow(d1) > 10)
})

test_that("Data outputing works for Plot_VarsByGroup", {
    expect_type(
        p <- dittoPlotVarsAcrossGroups(
            c("gene1", "gene2"), object = seurat, group.by = "age",
            data.out = TRUE),
        "list")
    expect_true(length(p) == 2)
    expect_type(p$data, "list")
})

test_that("Data outputing works for Heatmap", {
    expect_type(
        hm <- dittoHeatmap(c("gene1", "gene2"), object = seurat,
            data.out = TRUE),
        "list")
    expect_true(length(hm) == 2)
    expect_true("mat" %in% names(hm[[1]]))
})

