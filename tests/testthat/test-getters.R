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
    expect_type(meta("score", seurat, adjustment = "z-score"),
        "double")
    expect_equal(0, mean(
        meta("score", seurat, adjustment = "z-score")))
    expect_type(meta("score", seurat, adjustment = "relative.to.max"),
        "double")
    expect_equal(0:1, range(
        meta("score", seurat, adjustment = "relative.to.max")))
    expect_equal(factor(meta("score", seurat)),
        meta("score", seurat, adj.fxn = function(x) {factor(x)}))
})

test_that("metaLevels works for Seurat and SCE", {
    expect_type(groups <- metaLevels("groups", seurat),
        "character")
    expect_equal(groups, metaLevels("groups", sce))
})

test_that("meta and metaLevels give error when given a non-meta", {
    expect_error(meta("a", seurat),
        "is not a metadata of 'object'", fixed = TRUE)
    expect_error(metaLevels("a", seurat),
        "is not a metadata of 'object'", fixed = TRUE)
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
    expect_equal(0:1, range(
        gene("gene1", seurat, adjustment = "relative.to.max")))
    expect_equal(gene("gene1", seurat)+1,
        gene("gene1", seurat, adj.fxn = function(x) {x+1}))
})

test_that("gene gives error when given a non-meta", {
    expect_error(gene("a", seurat),
        "is not a gene of 'object'", fixed = TRUE)
})

test_that("gene, isGene, and getGenes work with swap.rownames",{
    
    swap_genes <- paste(rownames(sce), "symb", sep = "_")
    
    expect_equal(
        getGenes(sce, swap.rownames = "symbol"),
        swap_genes)
    
    expect_true(
        isGene("gene1_symb", sce, swap.rownames = "symbol"))
    
    expect_type(
        gene("gene1_symb", sce, swap.rownames = "symbol"),
        "double"
    )
    
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
        "is not a metadata or gene nor equal in length to ncol('object')", fixed = TRUE)
})

test_that("isBulk works properly", {
    expect_true(isBulk(bulk))
    expect_false(isBulk(sce))
    expect_false(isBulk(seurat))
    expect_true(isBulk(as(sce, "SummarizedExperiment")))
})

test_that("setBulk works properly", {
    expect_true(isBulk(setBulk(bulk)))
    expect_false(isBulk(setBulk(bulk, FALSE)))
    expect_true(isBulk(setBulk(sce)))
})

test_that(".which_cells converts non-string cells.use to string", {
    expect_equal(.which_cells(1:10, sce), colnames(sce)[1:10])
    logical <- rep(FALSE, ncol(sce))
    logical[1:10] <- TRUE
    expect_equal(.which_cells(logical, sce), colnames(sce)[1:10])
    expect_equal(.which_cells(colnames(sce)[1:10], sce), colnames(sce)[1:10])
})

test_that(".which_cells errors when logical 'cells.use' is the wrong length", {
    expect_error(.which_cells(TRUE, bulk),
        "'cells.use' length must equal the number of cells/samples in 'object' when given in logical form",
        fixed = TRUE)
})
