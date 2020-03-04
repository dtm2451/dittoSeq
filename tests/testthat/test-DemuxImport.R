# Tests for importDemux2Seurat function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DemuxImport.R")

test_that("importDemux2Seurat works when all barcodes are there", {
    expect_s4_class(
        importDemux2Seurat(
            seurat, lane.meta = "groups",
            Demuxlet.best = "mock_demux.best"),
        "Seurat")
})

test_that("importDemux2Seurat can name lanes with 'lane.names'", {
    expect_s4_class(
        t <- importDemux2Seurat(
            seurat, lane.meta = "groups",
            lane.names = meta.levels("groups", seurat),
            Demuxlet.best = "mock_demux.best"
            ),
        "Seurat")
    expect_equal(
        meta.levels("Lane",t),
        meta.levels("groups", seurat)
    )
})

test_that("importDemux2Seurat works when some barcodes are not there", {
    expect_s4_class(
        importDemux2Seurat(
            seurat, lane.meta = "groups", Demuxlet.best = "mock_demux_missing10.best"),
        "Seurat")
})

seurat2 <- seurat
colnames(seurat2@assays$RNA@data) <-
    paste(colnames(seurat2), rep(1:2, each = 40)[seq_len(ncells)], sep = "-")
test_that("importDemux2Seurat reads lanes from -# names", {
    expect_s4_class(
        t <- importDemux2Seurat(
            seurat2, Demuxlet.best = "mock_demux_-num2Lanes.best"),
        "Seurat")
    expect_equal(length(meta.levels("Lane",t)), 2)
})

test_that("importDemux2Seurat can be verbose or not", {
    expect_message(
        importDemux2Seurat(
            seurat, Demuxlet.best = "mock_demux.best",
            verbose = TRUE))
    expect_message(
        importDemux2Seurat(
            seurat, Demuxlet.best = "mock_demux.best",
            verbose = FALSE),
        NA)
})

seurat.demux <- importDemux2Seurat(
    seurat, lane.meta = "groups",
    Demuxlet.best = "mock_demux.best",
    verbose = FALSE)

test_that("demux.SNP.summary version of dittoPlot works", {
    expect_s3_class(
        demux.SNP.summary(seurat.demux),
        "ggplot")
})

test_that("demux.calls.summary works", {
    expect_s3_class(
        demux.calls.summary(seurat.demux),
        "ggplot")
})

test_that("demux.calls.summary options work", {
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            singlets.only = TRUE),
        "ggplot")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            rotate.labels = FALSE),
        "ggplot")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            data.out = TRUE),
        "data.frame")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            theme = theme_bw()),
        "ggplot")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            color = "yellow"),
        "ggplot")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            main = "1",
            sub = "2",
            xlab = "3",
            ylab = "4"),
        "ggplot")
    expect_s3_class(
        demux.calls.summary(seurat.demux,
            main = NULL,
            sub = NULL,
            xlab = NULL,
            ylab = NULL),
        "ggplot")
})
