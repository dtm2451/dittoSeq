# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Heatmap.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
seurat$number2 <- as.numeric(rev(seq_along(colnames(seurat))))
genes <- getGenes(seurat)[1:9]

test_that("Heatmap can be plotted for Seurat or SCE", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = sce),
        "pheatmap")
})

test_that("Heatmap title can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            main = "Hello there!"),
        "pheatmap")
})

test_that("Heatmap data type can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            slot = "counts"), # normally raw.data, but as.Seurat......
        "pheatmap")
})

test_that("Heatmap can hide rownames/colnames", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            show.colnames = FALSE),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            show.rownames = FALSE),
        "pheatmap")
})

test_that("Heatmap gives proper warnings when it should", {
    # Function throws a warning when any genes are not captured in the target cells.
    #   One of the genes are not expressed in the first 10 cells of seurat.
    expect_warning(
        dittoHeatmap(
            genes,
            object = seurat,
            cells.use = meta("number", seurat)<5),
        NULL)
    # Function throws an error when no genes provided are captured in the target cells.
    expect_error(
        dittoHeatmap(
            c("MS4A1","CD14","FCER1A"),
            object = seurat,
            cells.use = meta("number", seurat)<10),
        NULL)
})

test_that("Heatmap sample renaming by metadata works", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cell.names.meta = "number"),
        "pheatmap")
})

test_that("Heatmap highlight genes works", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            show.colnames = FALSE,
            show.rownames = FALSE,
            highlight.genes = "gene1"),
        "pheatmap")
})

test_that("Heatmap can be scaled to max", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            scaled.to.max = TRUE),
        "pheatmap")
})

test_that("Heatmap colors can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            heatmap.colors = colorRampPalette(c("yellow", "black", "red"))(50)),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            scaled.to.max = TRUE,
            heatmap.colors.max.scaled = colorRampPalette(c("black", "yellow"))(25)),
        "pheatmap")
})

test_that("Heatmap annotations can be given through metadata provision", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "ident"),
        "pheatmap")
})

test_that("Heatmap can be subset to certain cells", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cells.use = meta("number", seurat)<60),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cells.use = colnames(seurat)[meta("number", seurat)<60]),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "ident",
            cells.use = colnames(seurat)[meta("number", seurat)<60]),
        "pheatmap")
})

test_that("Heatmap annotation colors can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = c("number","ident"),
            annotation.colors = c("red", "yellow", "blue", "purple", "green3")),
        "pheatmap")
})

test_that("Coloring works for discrete column and row annotations", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "ident",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                genes,
                row.names = genes)),
        "pheatmap")
})

test_that("Coloring works for continuous column and row annotations", {
    # "Coloring works for discrete column annotations if BOTH PASS && LOOK DIFFERENT"
        # If: annotations are all continuous.
        # And...
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = genes),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
        # And if column annotation bar flips, but legend stays the same.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number2",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = genes),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
})

test_that("Heatmap can be ordered by metadata, expression, or user-input vector", {
    # Works for expression
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            order.by = "gene1"),
        "pheatmap")
    # Works for metadata
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            order.by = "groups",
            annotation.metas = "groups"),
        "pheatmap")
    # Works with vectors provided
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number",
            order.by = seq_along(colnames(seurat))),
        "pheatmap")
})

test_that("Heatmap can be ordered when also subset to certain cells", {
    # Works for expression
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            order.by = "gene1",
            cells.use = meta("number", seurat)<60),
        "pheatmap")
    # Works for metadata
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            order.by = "groups",
            annotation.metas = "groups",
            cells.use = colnames(seurat)[meta("number", seurat)<60]),
        "pheatmap")
    # Works with vectors provided
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number",
            order.by = seq_along(colnames(seurat)),
            cells.use = colnames(seurat)[meta("number", seurat)<60]),
        "pheatmap")
})
