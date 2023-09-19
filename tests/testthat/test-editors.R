# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); for (file in list.files("../../R", "utils", full.names = TRUE)) {source(file)}; source("test-editors.R")

sce_alts <- sce
ori <- rownames(sce_alts)
symbols <- paste0('symb.', rownames(sce_alts))
rowData(sce_alts)$symbol <- symbols
altExp(sce_alts) <- sce_alts[34:66,]
altExp(sce_alts, 'ADT') <- sce_alts[67:100,]
sce_alts <- sce_alts[1:33,]

test_that(".swap_rownames, explicit-method, updates rownames of main exp", {
    expect_equal(
        getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c(main = 'symbol')),
        c(symbols[1:33], ori[34:100])
    )
})

test_that(".swap_rownames, explicit-method, updates rownames of generic alternate exp", {
    expect_equal(
        getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c(altexp = 'symbol')),
        c(ori[1:33], symbols[34:66], ori[67:100])
    )
})

test_that(".swap_rownames, explicit-method, updates rownames of named alternate exp", {
    expect_equal(
        getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c(ADT = 'symbol')),
        c(ori[1:66], symbols[67:100])
    )
})

test_that(".swap_rownames, explicit-method, updates rownames of multiple targets", {
    expect_equal(
        getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c('symbol', altexp = 'symbol', ADT = 'symbol')),
        symbols
    )
})

test_that(".swap_rownames, simple-method, updates rownames of multiple targets", {

    expect_equal(
        getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'symbol'),
        symbols
    )

    sce_alts_mod1 <- sce_alts
    rowData(sce_alts_mod1)$other <- rowData(sce_alts_mod1)$symbol
    rowData(sce_alts_mod1)$symbol <- NULL

    sce_alts_mod2 <- sce_alts
    rowData(SingleCellExperiment::altExp(sce_alts_mod2))$other <- rowData(SingleCellExperiment::altExp(sce_alts_mod2))$symbol
    rowData(altExp(sce_alts_mod2))$symbol <- NULL

    # Also warn when no match for a main experiment
    expect_warning(
        getGenes(sce_alts_mod1, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'symbol'),
        "Skipping swapping rownames for main")
    expect_equal(
        suppressWarnings(getGenes(sce_alts_mod1, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'symbol')),
        c(ori[1:33], symbols[34:100])
    )
    # Secondary option used, no warnings because each has one
    expect_warning(
        getGenes(sce_alts_mod1, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c('symbol', 'other')),
        NA)
    expect_equal(
        getGenes(sce_alts_mod1, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c('symbol', 'other')),
        symbols
    )

    # Also warn when no match for an alternative experiment
    expect_warning(
        getGenes(sce_alts_mod2, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'symbol'),
        "Skipping swapping rownames for unnamed1")
    expect_equal(
        suppressWarnings(getGenes(sce_alts_mod2, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'symbol')),
        c(symbols[1:33], ori[34:66], symbols[67:100])
    )
    # Secondary option used, no warnings because each has one
    expect_warning(
        getGenes(sce_alts_mod2, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c('symbol', 'other')),
        NA)
    expect_equal(
        getGenes(sce_alts_mod2, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = c('symbol', 'other')),
        symbols
    )

    # No matches
    expect_equal(
        suppressWarnings(getGenes(sce_alts, assay = c('logcounts', 'altexp', 'ADT'), swap.rownames = 'nat')),
        ori
    )
})
