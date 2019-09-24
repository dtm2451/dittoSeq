# Tests for multi_dittoDimPlotVaryCells function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-VaryCells.R")

pbmc <- Seurat::pbmc_small
pbmc.se <- Seurat::as.SingleCellExperiment(pbmc)

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
grp <- "RNA_snn_res.0.8"
cont <- "CD14"
disc <- "RNA_snn_res.1"

test_that("VaryCells fxn can show continuous or discrete data", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp),
        "gtable")
    expect_s3_class(
        multi_dittoDimPlotVaryCells(disc, pbmc, grp),
        "gtable")
})

test_that("VaryCells fxn can output plots as a list", {
    expect_type(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            OUT.List = TRUE),
        "list")
})

test_that("VaryCells fxn can adjust how expression data is obtained", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            min = 0, max =2000),
        "gtable")
    #Manual Check: scales should be different in the next 2
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            data.type = "raw"),
        "gtable")
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp),
        "gtable")
})

test_that("VaryCells fxn errors as wanted when given 'cells.use'.", {
    expect_error(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            cells.use = colnames(pbmc)[1:5]),
        "Subsetting with cells.use is incopatible with this function.")
})

test_that("VaryCells fxn labels subsetting works", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            vary.cells.levels = 1:2),
        "gtable")
})

test_that("VaryCells 'show.' tweaks all work", {
    # Manual Check: Removes allcells and legend, and titles
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            show.allcells.plot = FALSE, show.legend.single = FALSE,
            show.titles = FALSE),
        "gtable")
    # Manual Check: Adds legends to all plots
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            show.legend.allcells.plot = TRUE, show.legend.plots = TRUE),
        "gtable")
})

test_that("VaryCells allcells title can be changed examples work", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, pbmc, grp,
            allcells.main = "DIFFERENT"),
        "gtable")
})

test_that("VaryCells color.panel can be adjusted", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(disc, pbmc, grp,
            color.panel = c("red","blue","yellow","gray50","purple"),
            colors = 5:3),
        "gtable")
})


