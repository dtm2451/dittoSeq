# Tests for multi_dittoDimPlotVaryCells function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-multi_VaryCells.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
grp <- "age"
cont <- "gene2"
disc <- "groups"

test_that("VaryCells fxn can show continuous or discrete data", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp),
        "gtable")
    expect_s3_class(
        multi_dittoDimPlotVaryCells(disc, object=seurat, grp),
        "gtable")
})

test_that("VaryCells works for SCE", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=sce, grp),
        "gtable")
})

test_that("VaryCells fxn can output plots as a list", {
    expect_type(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            list.out = TRUE),
        "list")
})

test_that("VaryCells fxn can adjust how expression data is obtained", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            min = 0, max =2000),
        "gtable")
    #Manual Check: scales should be different in the next 2
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            slot = "counts"),
        "gtable")
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp),
        "gtable")
})

test_that("VaryCells fxn levels subsetting works", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            vary.cells.levels = 1:2),
        "gtable")
})

test_that("VaryCells 'show.' tweaks all work", {
    # Manual Check: Removes allcells & legend & titles
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            show.allcells.plot = FALSE, show.legend.single = FALSE,
            show.titles = FALSE),
        "gtable")
    # Manual Check: Adds legends to all plots
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            show.legend.allcells.plot = TRUE, show.legend.plots = TRUE),
        "gtable")
})

test_that("VaryCells allcells title can be changed", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            allcells.main = "DIFFERENT"),
        "gtable")
})

test_that("VaryCells color.panel can be adjusted", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(disc, object=seurat, grp,
            color.panel = c("red","blue","yellow","gray50","purple"),
            colors = 5:1),
        "gtable")
})

test_that("VaryCells color.panel can be adjusted", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(disc, object=seurat, grp,
            color.panel = c("red","blue","yellow","gray50","purple"),
            colors = 5:1),
        "gtable")
})

test_that("VaryCells fxn errors as wanted when given 'cells.use'.", {
    expect_error(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            cells.use = colnames(seurat)[1:5]),
        "Further subsetting with 'cells.use'", fixed = TRUE)
})

test_that("VaryCells tells that 'main' is ignored.", {
    expect_message(
        multi_dittoDimPlotVaryCells(cont, object=seurat, grp,
            main = "HELLO"),
        "'main' ignored", fixed = TRUE)
})

test_that("VaryCells swap.rownames works", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(
            sce, "gene1_symb", grp, swap.rownames = "symbol"),
        "gtable")
})
