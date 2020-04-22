# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-data-output.R")

test_that("Data outputing works for ScatterPlot", {
    expect_type(
        data <- dittoScatterPlot("gene1", "gene2", object = seurat, data.out = TRUE),
        "list")
    expect_true(is.data.frame(data$Target_data))
    expect_true(is.data.frame(data$Others_data))
    expect_true(ncol(data$Target_data) >= 2 && nrow(data$Target_data) == ncol(seurat))
})

test_that("Data outputing works for DimPlot", {
    expect_type(
        data <- dittoDimPlot("gene1", object = seurat, data.out = TRUE),
        "list")
    expect_true(is.data.frame(data$Target_data))
    expect_true(is.data.frame(data$Others_data))
    expect_true(ncol(data$Target_data) > 2 && nrow(data$Target_data) == ncol(seurat))
})

test_that("Data outputing works for BarPlot", {
    expect_type(
        data <- dittoBarPlot("clusters", object = seurat, group.by = "age", data.out = TRUE),
        "list")
    expect_true(is.data.frame(data$data))
})

test_that("Data outputing works for Plot", {
    expect_type(
        data <- dittoPlot("gene1", object = seurat, group.by = "age", color.by = "age", data.out = TRUE),
        "list")
    expect_true(is.data.frame(data$data))
    expect_true(ncol(data$data) > 2 && nrow(data$data) == ncol(seurat))
})

test_that("Data outputing works for Plot_VarsByGroup", {
    expect_type(
        data <- dittoPlotVarsAcrossGroups(c("gene1", "gene2"), object = seurat, group.by = "age", data.out = TRUE),
        "list")
    expect_true(is.data.frame(data$data))
})

test_that("Data outputing works for Heatmap", {
    expect_type(
        hm <- dittoHeatmap(c("gene1", "gene2"), object = seurat,
            data.out = TRUE),
        "list")
    expect_true("mat" %in% names(hm))
})

