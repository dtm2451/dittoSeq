# Tests for visualization functions
# library(dittoSeq); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-getters.R")

test_that("getMetas works for Seurat.v3, SCE, and RNAseq", {
    expect_type(metas <- getMetas(pbmc),
        "character")
    expect_true(all(metas %in% getMetas(pbmc.rnaseq)))
    expect_true(all(metas %in% getMetas(pbmc.rnaseq)))
})

metas <- getMetas(pbmc)

test_that("isMeta works for Seurat.v3, SCE, and RNAseq", {
    expect_false(isMeta("HELLO", pbmc))
    expect_false(isMeta("HELLO", pbmc.se))
    expect_false(isMeta("HELLO", pbmc.rnaseq))
    expect_true(isMeta("ident", pbmc))
    expect_true(isMeta("nCount_RNA", pbmc))
    expect_true(isMeta("nCount_RNA", pbmc.se))
    expect_true(isMeta("nCount_RNA", pbmc.rnaseq))
})

test_that("meta works for Seurat.v3, SCE, and RNAseq", {
    expect_type(counts <- meta("nCount_RNA", pbmc),
        "double")
    expect_equal(counts, meta("nCount_RNA", pbmc.se))
    expect_equal(counts, meta("nCount_RNA", pbmc.rnaseq))
})

test_that("meta.levels works for Seurat.v3, SCE, and RNAseq", {
    expect_type(groups <- meta.levels("groups", pbmc),
        "character")
    expect_equal(groups, meta.levels("groups", pbmc.se))
    expect_equal(groups, meta.levels("groups", pbmc.rnaseq))
})

test_that("getGenes works for Seurat.v3, SCE, and RNAseq", {
    expect_type(genes <- getGenes(pbmc),
        "character")
    expect_equal(genes, getGenes(pbmc.se))
    expect_equal(genes, getGenes(pbmc.rnaseq))
})

genes <- getGenes(pbmc)

test_that("isGene works for Seurat.v3, SCE, and RNAseq", {
    expect_false(isGene("HELLO", pbmc))
    expect_false(isGene("HELLO", pbmc.se))
    expect_false(isGene("HELLO", pbmc.rnaseq))
    expect_true(isGene("MYL9", pbmc))
    expect_true(isGene("MYL9", pbmc.se))
    expect_true(isGene("MYL9", pbmc.rnaseq))
})

test_that("gene works for Seurat.v3, SCE, and RNAseq", {
    expect_type(MYL9 <- gene("MYL9", pbmc),
        "double")
    expect_equal(MYL9, gene("MYL9", pbmc.se))
    expect_equal(MYL9, gene("MYL9", pbmc.rnaseq))
})

test_that("gene works for different data.types (Seurat.v3)", {
    expect_type(gene("MYL9", pbmc, data.type = "raw"),
        "double")
    expect_type(gene("MYL9", pbmc, data.type = "normalized"),
        "double")
    expect_type(gene("MYL9", pbmc, data.type = "scaled"),
        "double")
    expect_type(gene("MYL9", pbmc, data.type = "relative"),
        "double")
    expect_equal(0, mean(
        gene("MYL9", pbmc, data.type = "relative")))
    expect_type(gene("MYL9", pbmc, data.type = "raw.normalized.to.max"),
        "double")
    expect_equal(1, max(
        gene("MYL9", pbmc, data.type = "raw.normalized.to.max")))
    expect_type(gene("MYL9", pbmc, data.type = "normalized.to.max"),
        "double")
    expect_equal(1, max(
        gene("MYL9", pbmc, data.type = "normalized.to.max")))
})

test_that("getReductions works for Seurat.v3, SCE, and RNAseq", {
    expect_type(reductions <- getReductions(pbmc),
        "character")
    expect_equal(reductions, tolower(getReductions(pbmc.se)))
    expect_equal(reductions, getReductions(pbmc.rnaseq))
})

