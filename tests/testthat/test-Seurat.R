# Tests for functionality with Seurat objects
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Seurat.R")

# Convert SCE to Seurat
seurat <- suppressWarnings(Seurat::as.Seurat(sce))
rownames(seurat@meta.data) <- colnames(seurat)

# For importDemux
sce.noDash <- sce
colnames(sce.noDash) <- sapply(colnames(sce.noDash), function(X) strsplit(X, "-")[[1]][1])

seurat.noDash <- suppressWarnings(Seurat::as.Seurat(sce.noDash))

sce.demux <- importDemux(
  sce.noDash, lane.meta = "groups",
  demuxlet.best = "mock_demux.best",
  verbose = FALSE)
seurat.demux <- suppressWarnings(Seurat::as.Seurat(sce.demux))

best <- read.table(file = "mock_demux.best", header=TRUE, sep="\t", stringsAsFactors = FALSE)
samples <- vapply(
  as.character(best$BEST),
  function(X) strsplit(X,'-')[[1]][2],
  FUN.VALUE = character(1))[seq_len(ncol(sce))]
names(samples) <- NULL

### TESTS
test_that("dittoBarPlot works for a Seurat", {
  expect_s3_class(
    dittoBarPlot(
      seurat, "clusters", group.by = "age"),
    "ggplot")
})

test_that("importDemux works for a Seurat", {
  expect_s4_class(
    t <- importDemux(
      object = seurat.noDash, lane.meta = "groups",
      demuxlet.best = "mock_demux.best"),
    "Seurat")
  # All called samples are correct
  expect_true(all(samples == meta("Sample", t)))
})

test_that("demux.SNP.summary & demux.calls.summary work for a Seurat", {
  expect_s3_class(
    demux.SNP.summary(seurat.demux),
    "ggplot")
  expect_s3_class(
    demux.calls.summary(seurat.demux),
    "ggplot")
})
  
test_that("dittoDimPlot work for a Seurat & 'slot' input usable", {
  expect_s3_class(
    dittoDimPlot(
      "groups", object=seurat),
    "ggplot")
  
  df <- dittoDimPlot("gene1", object = seurat, data.out = TRUE,
      slot = "counts")
  expect_equal(
    df$Target_data$color,
    round(df$Target_data$color,0))
})

test_that("dittoDotPlot works with gene and meta data for a Seurat, with 'slot' doing what it should", {
  expect_s3_class(
    print(dittoDotPlot(seurat, group.by = "clusters",
        getGenes(sce)[1:5])),
    "ggplot")
  expect_s3_class(
    dittoDotPlot(seurat, group.by = "clusters",
        c("score", "score2", "score3")),
    "ggplot")
  expect_s3_class(
    dittoDotPlot(seurat, group.by = "clusters",
        c("score", "gene1")),
    "ggplot")
  
  expect_type(
    d_raw <- dittoDotPlot(seurat, getGenes(sce)[1:5], "clusters", data.out = TRUE, scale = FALSE,
        slot = "counts"),
    "list")
  expect_type(
    d_log <- dittoDotPlot(seurat, getGenes(sce)[1:5], "clusters", data.out = TRUE, scale = FALSE,
        slot = "data"),
    "list")
  expect_true(all(
    d_raw$data$color >= d_log$data$color))
})

test_that("Heatmap can be plotted for a Seurat and slot adjusted", {
  expect_s3_class(
    dittoHeatmap(
        genes = getGenes(sce)[1:9],
        object = seurat),
    "pheatmap")
  expect_s3_class(
    dittoHeatmap(
      genes = getGenes(sce)[1:9],
      object = seurat,
      slot = "counts"),
    "pheatmap")
  expect_true(all(
    (dittoHeatmap(
      genes = getGenes(sce)[1:9], object = seurat, data.out = TRUE)$mat) <=
      (dittoHeatmap(
        genes = getGenes(sce)[1:9], object = seurat, data.out = TRUE, slot = "counts")$mat)
  ))
})