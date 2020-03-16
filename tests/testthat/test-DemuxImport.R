# Tests for importDemux function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DemuxImport.R")

test_that("importDemux works when all barcodes are there", {
    expect_s4_class(
        importDemux(
            seurat, lane.meta = "groups",
            demuxlet.best = "mock_demux.best"),
        "Seurat")
})

test_that("importDemux can name lanes with 'lane.names'", {
    expect_s4_class(
        t <- importDemux(
            seurat, lane.meta = "groups",
            lane.names = metaLevels("groups", seurat),
            demuxlet.best = "mock_demux.best"
            ),
        "Seurat")
    expect_equal(
        metaLevels("Lane",t),
        metaLevels("groups", seurat)
    )
})

test_that("importDemux works when some barcodes are not there", {
    expect_s4_class(
        importDemux(
            seurat, lane.meta = "groups", demuxlet.best = "mock_demux_missing10.best"),
        "Seurat")
})

seurat2 <- seurat
colnames(seurat2@assays$RNA@data) <-
    paste(colnames(seurat2), rep(1:2, each = 40)[seq_len(ncells)], sep = "-")
test_that("importDemux reads lanes from -# names", {
    expect_s4_class(
        t <- importDemux(
            seurat2, demuxlet.best = "mock_demux_-num2Lanes.best"),
        "Seurat")
    expect_equal(length(metaLevels("Lane",t)), 2)
})

test_that("importDemux can be verbose or not", {
    expect_message(
        importDemux(
            seurat, demuxlet.best = "mock_demux.best",
            verbose = TRUE))
    expect_message(
        importDemux(
            seurat, demuxlet.best = "mock_demux.best",
            verbose = FALSE),
        NA)
})

seurat.demux <- importDemux(
    seurat, lane.meta = "groups",
    demuxlet.best = "mock_demux.best",
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
