# Tests for importDemux function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DemuxImport.R")

sce.noDash <- sce
colnames(sce.noDash) <- sapply(colnames(sce.noDash), function(X) strsplit(X, "-")[[1]][1])
seurat.noDash <- suppressWarnings(Seurat::as.Seurat(sce.noDash))

sce2 <- sce.noDash
colnames(sce2) <- paste(colnames(sce2), rep(1:2, each = 40)[seq_len(ncells)], sep = "-")
seurat2 <- suppressWarnings(Seurat::as.Seurat(sce2))

sce3 <- sce.noDash
colnames(sce3) <- paste(colnames(sce3), c(rep(1, 45),rep(2, 35))[seq_len(ncells)], sep = "-")
seurat3 <- suppressWarnings(Seurat::as.Seurat(sce3))

### Extract the proper samples calls
best <- read.table(file = "mock_demux.best", header=TRUE, sep="\t", stringsAsFactors = FALSE)
samples <- vapply(
    as.character(best$BEST),
    function(X) strsplit(X,'-')[[1]][2],
    FUN.VALUE = character(1))[seq_len(ncol(sce))]
names(samples) <- NULL
# Vaariations for particular tests
samples_NAatEnd <- samples
samples_NAatEnd[71:ncol(sce)] <- NA

test_that("importDemux works when all barcodes are there (Seurat and SCE)", {
    expect_s4_class(
        t <- importDemux(
            object = seurat.noDash, lane.meta = "groups",
            demuxlet.best = "mock_demux.best"),
        "Seurat")
    expect_s4_class(
        importDemux(
            object = sce.noDash, lane.meta = "groups",
            demuxlet.best = "mock_demux.best"),
        "SingleCellExperiment")
    # All called samples are correct
    expect_true(all(samples == meta("Sample", t)))
})

test_that("importDemux can name lanes with 'lane.names' (Seurat and SCE)", {
    expect_s4_class(
        t <- importDemux(
            seurat.noDash, lane.meta = "groups",
            lane.names = metaLevels("groups", seurat),
            demuxlet.best = "mock_demux.best"
            ),
        "Seurat")
    expect_equal(
        metaLevels("Lane",t),
        metaLevels("groups", seurat)
    )
    expect_s4_class(
        t <- importDemux(
            sce.noDash, lane.meta = "groups",
            lane.names = metaLevels("groups", seurat),
            demuxlet.best = "mock_demux.best"
            ),
        "SingleCellExperiment")
    # All called samples are correct
    expect_true(all(samples == meta("Sample", t)))
})

test_that("importDemux works when some barcodes are not there (Seurat and SCE)", {
    expect_s4_class(
        importDemux(
            seurat.noDash, lane.meta = "groups", demuxlet.best = "mock_demux_missing10.best"),
        "Seurat")
    expect_s4_class(
        t <- importDemux(
            sce.noDash, lane.meta = "groups", demuxlet.best = "mock_demux_missing10.best"),
        "SingleCellExperiment")
    # All called samples are correct
    expect_true(all(samples_NAatEnd == meta("Sample", t), na.rm = TRUE))
    # Non-called samples are the missing cells
    expect_true(all(is.na(samples_NAatEnd) == is.na(meta("Sample", t))))
})

test_that("importDemux works with incremented -# and multiple .best given (Seurat and SCE)", {
    expect_s4_class(
        importDemux(
            seurat3,
            demuxlet.best = c("mock_demux_1st45.best",
                              "mock_demux_2nd35.best")),
        "Seurat")
    expect_s4_class(
        t <- importDemux(
            sce3,
            demuxlet.best = c("mock_demux_1st45.best",
                              "mock_demux_2nd35.best")),
        "SingleCellExperiment")
    # All called samples are correct
    expect_true(all(samples == meta("Sample", t)))
})

test_that("importDemux works with incremented -# and multiple .best given when some cells not in a .best (Seurat and SCE)", {
    expect_s4_class(
        importDemux(
            seurat2,
            demuxlet.best = c("mock_demux_1st45.best",
                              "mock_demux_2nd35.best")),
        "Seurat")
    expect_s4_class(
        t <- importDemux(
            sce2,
            demuxlet.best = c("mock_demux_1st45.best",
                              "mock_demux_2nd35.best")),
        "SingleCellExperiment")
    # All called samples are correct
    expect_true(all(samples[!is.na(meta("Sample", t))] == meta("Sample", t)[!is.na(meta("Sample", t))]))
    # Non-called samples are the 5 with non-matching lanes
    expect_true(all(is.na(meta("Sample", t)[41:45])))
})

test_that("importDemux reads lanes from -# names", {
    expect_s4_class(
        t <- importDemux(
            seurat2, demuxlet.best = "mock_demux_-num2Lanes.best"),
        "Seurat")
    expect_equal(length(metaLevels("Lane",t)), 2)
    expect_s4_class(
        t <- importDemux(
            sce2, demuxlet.best = "mock_demux_-num2Lanes.best"),
        "SingleCellExperiment")
    expect_equal(length(metaLevels("Lane",t)), 2)
    # All called samples are correct
    expect_true(all(samples_NAatEnd == meta("Sample", t), na.rm = TRUE))
    # Non-called samples are the expected missing cells
    expect_true(all(is.na(samples_NAatEnd) == is.na(meta("Sample", t))))
})

test_that("importDemux can be verbose or not", {
    expect_message(
        importDemux(
            seurat.noDash, demuxlet.best = "mock_demux.best",
            verbose = TRUE))
    expect_message(
        importDemux(
            seurat.noDash, demuxlet.best = "mock_demux.best",
            verbose = FALSE),
        NA)
})

test_that("importDemux gives error when no barcodes match", {
    expect_error(
        importDemux(
            seurat, demuxlet.best = "mock_demux.best",
            verbose = FALSE),
        "No barcodes match between 'object' and 'demuxlet.best'", fixed = TRUE)
})

seurat.demux <- importDemux(
    seurat.noDash, lane.meta = "groups",
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
