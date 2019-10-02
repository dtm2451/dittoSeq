# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-hover.R")

test_that("Showing hover.data works for ScatterPlot", {
    expect_s3_class(
        dittoScatterPlot("MS4A1", "GNLY", object = pbmc, do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident")),
        "plotly")
})

test_that("Showing hover.data works for DimPlot", {
    expect_s3_class(
        dittoDimPlot("MS4A1", object = pbmc, do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident")),
        "plotly")
    ### Manual Check: gene counts should become integers
    expect_s3_class(
        dittoDimPlot("MS4A1", object = pbmc, do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident"),
            hover.data.type = "raw"),
        "plotly")
})

test_that("Showing hover.data works for BarPlot", {
    expect_s3_class(
        dittoBarPlot(
            "RNA_snn_res.0.8", object = pbmc,
            group.by = "RNA_snn_res.1",
            do.hover = TRUE),
        "plotly")
})

test_that("Showing hover.data works for Plot", {
    expect_s3_class(
        dittoPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident")),
        "plotly")
    expect_s3_class(
        dittoBoxPlot(
            "MS4A1", object = pbmc,
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident")),
        "plotly")
    expect_warning(
        dittoRidgePlot(
            "MS4A1", object = pbmc, plots = c("ridgeplot", "jitter"),
            group.by = "RNA_snn_res.1", color.by = "RNA_snn_res.1",
            do.hover = TRUE,
            hover.data = c("MS4A1","RNA_snn_res.0.8","ident")),
        NULL)
})

test_that("Showing hover.data works for VarsAcrossGroups", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(c("MS4A1", "GNLY", "CD3E"),
            object = pbmc, group.by = "RNA_snn_res.1",
            do.hover = TRUE),
        "plotly")
})
