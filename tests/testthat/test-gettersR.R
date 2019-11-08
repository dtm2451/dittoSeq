# Tests for visualization functions
# library(dittoSeq); library(testthat); source("tests/testthat/setup.R"); source("tests/testthat/test-getters.R")

test_that("get.metas works for Seurat.v3, SCE, and RNAseq", {
    expect_type(metas <- get.metas(pbmc),
        "character")
    expect_true(all(metas %in% get.metas(pbmc.rnaseq)))
    expect_true(all(metas %in% get.metas(pbmc.rnaseq)))
})

metas <- get.metas(pbmc)

test_that("is.meta works for Seurat.v3, SCE, and RNAseq", {
    expect_false(is.meta("HELLO", pbmc))
    expect_false(is.meta("HELLO", pbmc.se))
    expect_false(is.meta("HELLO", pbmc.rnaseq))
    expect_true(is.meta("ident", pbmc))
    expect_true(is.meta("nCount_RNA", pbmc))
    expect_true(is.meta("nCount_RNA", pbmc.se))
    expect_true(is.meta("nCount_RNA", pbmc.rnaseq))
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

test_that("get.genes works for Seurat.v3, SCE, and RNAseq", {
    expect_type(genes <- get.genes(pbmc),
        "character")
    expect_equal(genes, get.genes(pbmc.se))
    expect_equal(genes, get.genes(pbmc.rnaseq))
})

genes <- get.genes(pbmc)

test_that("is.gene works for Seurat.v3, SCE, and RNAseq", {
    expect_false(is.gene("HELLO", pbmc))
    expect_false(is.gene("HELLO", pbmc.se))
    expect_false(is.gene("HELLO", pbmc.rnaseq))
    expect_true(is.gene("MYL9", pbmc))
    expect_true(is.gene("MYL9", pbmc.se))
    expect_true(is.gene("MYL9", pbmc.rnaseq))
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

test_that("get.reductions works for Seurat.v3, SCE, and RNAseq", {
    expect_type(reductions <- get.reductions(pbmc),
        "character")
    expect_equal(reductions, tolower(get.reductions(pbmc.se)))
    expect_equal(reductions, get.reductions(pbmc.rnaseq))
})

