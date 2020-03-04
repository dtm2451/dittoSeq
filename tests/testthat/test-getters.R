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

test_that("metaLevels works for Seurat and SCE", {
    expect_type(groups <- metaLevels("groups", seurat),
        "character")
    expect_equal(groups, metaLevels("groups", sce))
})

test_that("meta and metaLevels give error when given a non-meta", {
    expect_error(meta("a", seurat),
        "\"a\" is not a metadata of 'object'", fixed = TRUE)
    expect_error(metaLevels("a", seurat),
        "\"a\" is not a metadata of 'object'", fixed = TRUE)
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

test_that("gene gives error when given a non-meta", {
    expect_error(gene("a", seurat),
        "\"a\" is not a gene of 'object'", fixed = TRUE)
})

test_that("getReductions works for Seurat and SCE", {
    expect_type(reductions <- getReductions(seurat),
        "character")
    expect_equal(reductions, getReductions(sce))
})

test_that(".var_or_get_meta_or_gene gets metas, genes, spits back var, or errors if wrong length", {
    expect_equal(
        .var_OR_get_meta_or_gene("groups", seurat),
        meta("groups", seurat))
    expect_equal(
        .var_OR_get_meta_or_gene("gene1", seurat),
        gene("gene1", seurat))
    expect_equal(
        length(.var_OR_get_meta_or_gene(seq_len(ncol(seurat)), seurat)),
        ncol(seurat)) # just checks length because names are added by the function.
    expect_error(.var_OR_get_meta_or_gene(1,sce),
        "'var' is not a metadata or gene nor equal in length to ncol('object')", fixed = TRUE)
})


