# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); for (file in list.files("../../R", "utils", full.names = TRUE)) {source(file)}; source("../../R/gene-getters.R"); source("test-getters-multi-assay.R")

sce_alts <- sce
ori <- rownames(sce_alts)
symbols <- paste0('symb.', rownames(sce_alts))
rowData(sce_alts)$symbol <- symbols
altExp(sce_alts) <- sce_alts[34:66,]
altExp(sce_alts, 'ADT') <- sce_alts[67:100,]
sce_alts <- sce_alts[1:33,]

# Make Seurat, if can
try({
    seurat <- Seurat::as.Seurat(sce[1:33,])
    suppressWarnings(seurat[['RNA']] <- Seurat::as.Seurat(sce[34:66,])[['originalexp']])
    suppressWarnings(seurat[['ADT']] <- Seurat::as.Seurat(sce[67:100,])[['originalexp']])
    }, silent = TRUE)
seurat_conversion_worked <- exists("seurat")

# Ensure metadata has row.names (a past issue with the fxn)
if (seurat_conversion_worked) {
    rownames(seurat@meta.data) <- colnames(seurat)
}

test_that(".which_data works for all SCE alternate experiment direction methods", {
    expect_type(simple <- .which_data(c('logcounts', 'altexp', 'ADT'), NA, sce_alts),
                "double")
    expect_type(explicit <- .which_data(c('logcounts', 'altexp' = 'counts', 'ADT'), NA, sce_alts),
                "double")
    expect_equal(simple, explicit)
})

test_that("getGenes works across assays for Seurat and altExps for SCE", {
    expect_type(genes <- getGenes(sce_alts, assay = c('logcounts', 'altexp' = 'logcounts', 'ADT')),
                "character")
    expect_equal(genes, ori)

    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_equal(genes, getGenes(seurat, assay = c('originalexp', 'RNA', 'ADT')))
})

test_that("getGenes works across assays for Seurat and altExps for SCE", {
    expect_type(genes <- getGenes(sce_alts, assay = c('logcounts', 'altexp' = 'logcounts', 'ADT')),
        "character")
    expect_equal(genes, ori)

    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_equal(genes, getGenes(seurat, assay = c('originalexp', 'RNA', 'ADT')))
})

test_that("isGene works across assays for Seurat and altExps for SCE", {
    expect_equal(
        c(TRUE,FALSE),
        isGene(c("gene1", "gene100"), sce_alts, assay = c('logcounts', 'altexp')))
    expect_true(all(isGene("gene100", sce_alts, assay = c('logcounts', 'altexp', 'ADT'))))

    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_equal(
        c(TRUE,FALSE),
        isGene(c("gene1", "gene100"), seurat, assay = c('originalexp', 'RNA')))
    expect_true(all(isGene(c("gene1", "gene100"), seurat, assay = c('originalexp', 'RNA', 'ADT'))))
})

test_that("assay specification for SCE altExps works as intended", {
    # 'altexp' generalizes to first altExp
    expect_type(raw <- gene("gene50", sce_alts, assay = c('counts', 'altexp', ADT = 'logcounts')),
                "double")
    first <- SingleCellExperiment::altExpNames(sce_alts)[1]
    expect_type(raw2 <- gene("gene50", sce_alts, assay = first),
                "double")
    expect_equal(raw, raw2)

    # when only given altExp name, grabs from first slot
    expect_equal(
        as.vector(assay(altExp(sce_alts, first), 'counts')['gene50',]),
        as.vector(raw)
    )

    # when given as c(altexp = assay), grabs from chosen assay of the altExp
    expect_type(log <- gene("gene50", sce_alts, assay = c('altexp' = 'logcounts')),
                "double")
    expect_equal(
        as.vector(assay(altExp(sce_alts, first), 'logcounts')['gene50',]),
        as.vector(log)
    )
})

test_that(".var_or_get_meta_or_gene works across assays for Seurat and altExps for SCE", {
    expect_equal(
        .var_OR_get_meta_or_gene("groups", sce_alts, assay = c('logcounts', 'altexp' = 'logcounts')),
        meta("groups", sce_alts))
    expect_equal(
        .var_OR_get_meta_or_gene("gene50", sce_alts, assay = c('logcounts', 'altexp' = 'logcounts')),
        gene("gene50", sce_alts, assay = c('logcounts', 'altexp' = 'logcounts')))
    expect_equal(
        seq_len(ncol(sce_alts)),
        as.vector(.var_OR_get_meta_or_gene(
            seq_len(ncol(sce_alts)), sce_alts, assay = c('logcounts', 'altexp' = 'logcounts'))))

    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_equal(
        .var_OR_get_meta_or_gene("groups", seurat, assay = c('originalexp', 'RNA', 'ADT')),
        meta("groups", seurat))
    expect_equal(
        .var_OR_get_meta_or_gene("gene50", seurat, assay = c('originalexp', 'RNA', 'ADT')),
        gene("gene50", seurat, assay = c('originalexp', 'RNA', 'ADT')))
    expect_equal(
        seq_len(ncol(seurat)),
        as.vector(.var_OR_get_meta_or_gene(
            seq_len(ncol(seurat)), seurat, assay = c('originalexp', 'RNA', 'ADT'))))
})
