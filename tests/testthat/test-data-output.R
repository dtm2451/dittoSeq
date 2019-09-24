# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-data-output.R")

test_that("Data outputing works for ScatterPlot", {
    expect_type(
        dittoScatterPlot("MS4A1", "GNLY", object = pbmc,
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoScatterPlot("MS4A1", "GNLY", object = pbmc,
            data.out = TRUE)$Target_data,
        "list")
    expect_true(ncol(d1) >= 2 && nrow(d1) > 10)
})

test_that("Data outputing works for DimPlot", {
    expect_type(
        dittoDimPlot("MS4A1", object = pbmc,
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoDimPlot("MS4A1", object = pbmc,
            data.out = TRUE)$Target_data,
        "list")
    expect_true(ncol(d1) > 2 && nrow(d1) > 10)
})

test_that("Data outputing works for BarPlot", {
    expect_s3_class(
        dittoBarPlot(
            "RNA_snn_res.0.8", object = pbmc,
            group.by = "RNA_snn_res.1", data.out = TRUE,
            do.hover = TRUE),
        "data.frame")
})

test_that("Data outputing works for Plot", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_type(
        dittoPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE),
        "list")
    expect_type(
        d1 <- dittoPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE)$data,
        "list")
    expect_type(
        dittoRidgePlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE),
        "list")
    expect_type(
        d2 <- dittoRidgePlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE)$data,
        "list")
    expect_type(
        dittoBoxPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE),
        "list")
    expect_type(
        d3 <- dittoBoxPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            data.out = TRUE)$data,
        "list")
    expect_true(ncol(d1) > 2 && nrow(d1) > 10)
    expect_true(ncol(d2) > 2 && nrow(d1) > 10)
    expect_true(ncol(d3) > 2 && nrow(d1) > 10)
})

test_that("Data outputing works for Plot_VarsByGroup", {
    expect_type(
        p <- dittoPlotVarsAcrossGroups(
            c("MS4A1", "GNLY"), pbmc, group.by = "RNA_snn_res.1",
            data.out = TRUE),
        "list")
    expect_true(length(p) == 2)
    expect_type(p$data, "list")
})

test_that("Data outputing works for Heatmap", {
    expect_type(
        hm <- dittoHeatmap(c("MS4A1", "GNLY"), pbmc,
            data.out = TRUE),
        "list")
    expect_true(length(hm) == 2)
    expect_true("mat" %in% names(hm[[1]]))
})

