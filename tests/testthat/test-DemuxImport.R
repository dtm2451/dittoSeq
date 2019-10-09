# Tests for importDemux2Seurat function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DemuxImport.R")

test_that("importDemux2Seurat works when all barcodes are there", {
    expect_s4_class(
        importDemux2Seurat(
            pbmc, lane.meta = "groups",
            Demuxlet.best = "mock_demux.best"),
        "Seurat")
})

test_that("importDemux2Seurat can name lanes with 'lane.names'", {
    expect_s4_class(
        t <- importDemux2Seurat(
            pbmc, lane.meta = "groups", lane.names = c("test1","test2"),
            Demuxlet.best = "mock_demux.best"
            ),
        "Seurat")
    expect_equal(
        meta.levels("Lane",t),
        c("test1","test2")
    )
})

test_that("importDemux2Seurat works when some barcodes are not there", {
    expect_s4_class(
        importDemux2Seurat(
            pbmc, lane.meta = "groups", Demuxlet.best = "mock_demux_missing10.best"),
        "Seurat")
})

colnames(pbmc@assays$RNA@data) <- paste(colnames(pbmc), rep(1:2, each = 40), sep = "-")
test_that("importDemux2Seurat reads lanes from -# names", {
    expect_s4_class(
        t <- importDemux2Seurat(
            pbmc, Demuxlet.best = "mock_demux_-num2Lanes.best"),
        "Seurat")
    expect_equal(length(meta.levels("Lane",t)), 2)
})
