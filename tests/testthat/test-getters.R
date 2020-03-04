# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-getters.R")

test_that("getMetas works for Seurat and SCE", {
    expect_type(metas <- getMetas(seurat),
        "character")
    expect_true(all(metas %in% getMetas(sce)))
})

metas <- getMetas(seurat)

test_that("isMeta works for Seurat and SCE", {
    expect_false(isMeta("HELLO", seurat))
    expect_false(isMeta("HELLO", sce))
    expect_true(isMeta("ident", seurat))
    expect_true(isMeta("score", seurat))
    expect_true(isMeta("score", sce))
})

test_that("meta works for Seurat and SCE", {
    expect_type(counts <- meta("score", seurat),
        "double")
    expect_equal(counts, meta("score", sce))
})

test_that("meta.levels works for Seurat and SCE", {
    expect_type(groups <- meta.levels("groups", seurat),
        "character")
    expect_equal(groups, meta.levels("groups", sce))
})

test_that("getGenes works for Seurat and SCE", {
    expect_type(genes <- getGenes(seurat),
        "character")
    expect_equal(genes, getGenes(sce))
})

genes <- getGenes(seurat)

test_that("isGene works for Seurat and SCE", {
    expect_false(isGene("HELLO", seurat))
    expect_false(isGene("HELLO", sce))
    expect_true(isGene("gene1", seurat))
    expect_true(isGene("gene1", sce))
})

test_that("isGene works for different assays in SCE", {
    expect_false(isGene("HELLO", sce, assay = "counts"))
    expect_true(isGene("gene1", sce, assay = "counts"))
})

test_that("gene works for different data types for Seurat and SCE", {
    expect_type(gene("gene1", sce),
        "double")
    expect_type(gene("gene1", sce, assay = "counts"),
        "double")
    expect_type(gene("gene1", seurat, slot = "counts"),
        "double")
    expect_type(gene("gene1", seurat, slot = "data"),
        "double")
    expect_type(gene("gene1", seurat, adjustment = "z-score"),
        "double")
    expect_equal(0, mean(
        gene("gene1", seurat, adjustment = "z-score")))
    expect_type(gene("gene1", seurat, adjustment = "relative.to.max"),
        "double")
    expect_equal(1, max(
        gene("gene1", seurat, adjustment = "relative.to.max")))
})

test_that("getReductions works for Seurat and SCE", {
    expect_type(reductions <- getReductions(seurat),
        "character")
    expect_equal(reductions, getReductions(sce))
})

