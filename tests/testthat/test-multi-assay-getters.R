# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); for (file in list.files("../../R", "utils", full.names = TRUE)) {source(file)}; source("../../R/gene-getters.R"); source("test-multi-assay-getters.R")

sce_alts <- sce
ori <- rownames(sce_alts)
symbols <- paste0('symb.', rownames(sce_alts))
rowData(sce_alts)$symbol <- symbols
altExp(sce_alts) <- sce_alts[34:66,]
altExp(sce_alts, 'ADT') <- sce_alts[67:100,]
sce_alts <- sce_alts[1:33,]
# And a rowData that's only in the primary assay
rowData(sce_alts)$id <- paste0('id.', rownames(sce_alts))

# Additional altExp where:
# - new altexp name overlaps with an assay name (logcounts)
# - data are 0s or 1s for easy confirmation of grabbed data
sce_alts_overlaps <- sce_alts
new_alt <- sce_alts[1:33,]
counts(new_alt, withDimnames=FALSE) <- matrix(data = 0, nrow = nrow(new_alt), ncol = ncol(new_alt))
logcounts(new_alt, withDimnames=FALSE) <- matrix(data = 1, nrow = nrow(new_alt), ncol = ncol(new_alt))
altExp(sce_alts_overlaps, 'logcounts') <- new_alt

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

test_that("multi-assay edge cases", {
    # 'main' optionally specifies top-level data
    expect_true(
        identical(
            gene("gene1", sce_alts, assay = c(main='logcounts')),
            gene("gene1", sce_alts, assay = c('logcounts'))
        )
    )
    # 'altexp' alone specifies first-altExp's first assay
    expect_true(
        identical(
            gene("gene34", sce_alts, assay = c('altexp')),
            gene("gene34", sce_alts, assay = c('unnamed1'='counts'))
        )
    )

    # Priority to primary modality assay if overlap with an altExp name, & only grabs one
    expect_true(
        identical(
            gene("gene1", sce_alts_overlaps, assay = 'logcounts'),
            gene("gene1", sce_alts_overlaps, assay = c(main='logcounts'))
        )
    )
    expect_true(
        identical(
            gene("gene1", sce_alts_overlaps, assay = c('logcounts', logcounts='counts')),
            gene("gene1", sce_alts_overlaps, assay = c(main='logcounts'))
        )
    )
    expect_false(
        identical(
            gene("gene1", sce_alts_overlaps, assay = c('logcounts', logcounts='counts')),
            gene("gene1", sce_alts_overlaps, assay = c(logcounts='counts'))
        )
    )
})

test_that("multi-assay swap.rownames works as expected", {
    # Simple method: swaps for all modalities with the given rowData
    expect_true(
        all(
            isGene(c("symb.gene1", "symb.gene50"), object = sce_alts, assay = c('main', 'altexp'), swap.rownames = "symbol")
        )
    )
    # Simple method: swaps for any modality with the given rowData, warning for modalities without any
    expect_true(
        all(suppressWarnings(
            isGene(c("id.gene1", "gene50"), object = sce_alts, assay = c('main', 'altexp'), swap.rownames = "id")
        ))
    )
    # Simple method: prioritizes in given order
    expect_true(
        all(
            isGene(c("id.gene1", "symb.gene50"), object = sce_alts, assay = c('main', 'altexp'), swap.rownames = c("id", "symbol"))
        )
    )
    # Simple method warns when no swap.rownames found
    expect_warning(
        isGene(c("id.gene1", "gene50"), object = sce_alts, assay = c('main', 'altexp'), swap.rownames = "id"),
        "Skipping swapping rownames for unnamed1"
    )

    # Explicit method
    expect_true(
        all(
            isGene(c("id.gene1", "symb.gene50", "gene100"), object = sce_alts, assay = c('main', 'altexp', 'ADT'), swap.rownames = c(main='id', altexp="symbol"))
        )
    )
})

test_that("dittoPlot can make use of multi-assay feature access", {
    expect_s3_class(
        dittoPlot(
            c("id.gene1", "symb.gene50", "gene100"),
            object = sce_alts, group.by = "groups",
            assay = c('main', 'altexp', 'ADT'),
            swap.rownames = c(main='id', altexp="symbol")
        ),
        "ggplot")
    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_s3_class(
        dittoPlot(
            c("gene1", "gene50", "gene100"),
            object = seurat, group.by = "groups",
            assay = c('originalexp', 'RNA', 'ADT')
        ),
        "ggplot")
})

test_that("dittoDimPlot can make use of multi-assay feature access", {
    expect_s3_class(
        dittoDimPlot(
            c("id.gene1", "symb.gene50", "gene100"),
            object = sce_alts,
            assay = c('main', 'altexp', 'ADT'),
            swap.rownames = c(main='id', altexp="symbol")
        ),
        "ggplot")
    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_s3_class(
        dittoDimPlot(
            c("gene1", "gene50", "gene100"),
            object = seurat,
            assay = c('originalexp', 'RNA', 'ADT')
        ),
        "ggplot")
})

test_that("dittoDotPlot can make use of multi-assay feature access", {
    expect_s3_class(
        dittoDotPlot(
            c("id.gene1", "symb.gene50", "gene100"),
            object = sce_alts, group.by = "groups",
            assay = c('main', 'altexp', 'ADT'),
            swap.rownames = c(main='id', altexp="symbol")
        ),
        "ggplot")
    skip_if_not(seurat_conversion_worked, message = "Seurat conversion bug")
    expect_s3_class(
        dittoDotPlot(
            c("gene1", "gene50", "gene100"),
            object = seurat, group.by = "groups",
            assay = c('originalexp', 'RNA', 'ADT')
        ),
        "ggplot")
})

test_that("multi-assay documentation simple-path examples not confirmed elsewhere", {
    # assay = c('main', 'altexp', 'hto')
    expect_true(
        identical(
            .which_data(object = sce_alts, assay = c('main', 'altexp', 'ADT')),
            .which_data(object = sce_alts, assay = c(main='counts', altexp='counts', ADT='counts'))
        )
    )
    # assay = c('logexp', 'adt'='clr')
    expect_true(
        identical(
            .which_data(object = sce_alts, assay = c('logcounts', 'ADT'='logcounts')),
            .which_data(object = sce_alts, assay = c('main'='logcounts', ADT='logcounts'))
        )
    )
})
