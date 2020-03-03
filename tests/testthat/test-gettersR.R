# Tests for visualization functions
# library(dittoSeq); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-getters.R")

test_that("getMetas works for Seurat and SCE", {
    expect_type(metas <- getMetas(pbmc),
        "character")
    expect_true(all(metas %in% getMetas(pbmc.se)))
})

metas <- getMetas(pbmc)

test_that("isMeta works for Seurat and SCE", {
    expect_false(isMeta("HELLO", pbmc))
    expect_false(isMeta("HELLO", pbmc.se))
    expect_true(isMeta("ident", pbmc))
    expect_true(isMeta("nCount_RNA", pbmc))
    expect_true(isMeta("nCount_RNA", pbmc.se))
})

test_that("meta works for Seurat and SCE", {
    expect_type(counts <- meta("nCount_RNA", pbmc),
        "double")
    expect_equal(counts, meta("nCount_RNA", pbmc.se))
})

test_that("meta.levels works for Seurat and SCE", {
    expect_type(groups <- meta.levels("groups", pbmc),
        "character")
    expect_equal(groups, meta.levels("groups", pbmc.se))
})

test_that("getGenes works for Seurat and SCE", {
    expect_type(genes <- getGenes(pbmc),
        "character")
    expect_equal(genes, getGenes(pbmc.se))
})

genes <- getGenes(pbmc)

test_that("isGene works for Seurat and SCE", {
    expect_false(isGene("HELLO", pbmc))
    expect_false(isGene("HELLO", pbmc.se))
    expect_true(isGene("MYL9", pbmc))
    expect_true(isGene("MYL9", pbmc.se))
})

test_that("isGene works for different assays in SCE", {
    expect_false(isGene("HELLO", pbmc.se, assay = "counts"))
    expect_true(isGene("MYL9", pbmc.se, assay = "counts"))
})

test_that("gene works for different data types for Seurat and SCE", {
    expect_type(gene("MYL9", pbmc.se),
        "double")
    expect_type(gene("MYL9", pbmc.se, assay = "counts"),
        "double")
    expect_type(gene("MYL9", pbmc, slot = "counts"),
        "double")
    expect_type(gene("MYL9", pbmc, slot = "data"),
        "double")
    expect_type(gene("MYL9", pbmc, slot = "scale.data"),
        "double")
    expect_type(gene("MYL9", pbmc, adjustment = "z-score"),
        "double")
    expect_equal(0, mean(
        gene("MYL9", pbmc, adjustment = "z-score")))
    expect_type(gene("MYL9", pbmc, adjustment = "relative.to.max"),
        "double")
    expect_equal(1, max(
        gene("MYL9", pbmc, adjustment = "relative.to.max")))
})

test_that("getReductions works for Seurat and SCE", {
    expect_type(reductions <- getReductions(pbmc),
        "character")
    expect_equal(reductions, tolower(getReductions(pbmc.se)))
})

