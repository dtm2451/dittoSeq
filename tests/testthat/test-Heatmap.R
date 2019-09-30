# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Heatmap.R")

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
pbmc@meta.data$number2 <- as.numeric(rev(seq_along(colnames(pbmc))))
genes <- c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A","LYZ","PPBP","CD8A")

test_that("Heatmap can be plotted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc),
        "pheatmap")
})

test_that("Heatmap title can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            main = "Hello there!"),
        "pheatmap")
})

test_that("Heatmap data type can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            data.type = "raw"),
        "pheatmap")
})

test_that("Heatmap can hide rownames/colnames", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            show.colnames = FALSE),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            show.rownames = FALSE),
        "pheatmap")
})

test_that("Heatmap can be subset to certain cells", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            cells.use = meta("number", pbmc)<60),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            cells.use = colnames(pbmc)[meta("number", pbmc)<60]),
        "pheatmap")
})

test_that("Heatmap gives proper warnings when it should", {
    # Function gives and error if no geens are provided.
    expect_error(
        dittoHeatmap(
            c(),
            object = pbmc,
            cells.use = meta("number", pbmc)<10),
        NULL)
    # Function throws a warning when any genes are not captured in the target cells.
    #   Three of the genes are not expressed in the first 10 cells of pbmc.
    expect_warning(
        dittoHeatmap(
            genes,
            object = pbmc,
            cells.use = meta("number", pbmc)<10),
        NULL)
    # Function throws an error when no genes provided are captured in the target cells.
    expect_error(
        dittoHeatmap(
            c("MS4A1","CD14","FCER1A"),
            object = pbmc,
            cells.use = meta("number", pbmc)<10),
        NULL)
})

test_that("Heatmap highlight genes works", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            highlight.genes = "CD3E"),
        "pheatmap")
})

test_that("Heatmap can be scaled to max", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            scaled.to.max = TRUE),
        "pheatmap")
})

test_that("Heatmap colors can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            heatmap.colors = colorRampPalette(c("yellow", "black", "red"))(50)),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            scaled.to.max = TRUE,
            heatmap.colors.max.scaled = colorRampPalette(c("black", "yellow"))(25)),
        "pheatmap")
})

test_that("Heatmap annotations can be given through metadata provision", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
            annotation.metas = "ident"),
        "pheatmap")
})

test_that("Heatmap annotation colors can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = pbmc,
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
            object = pbmc,
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
            object = pbmc,
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
            object = pbmc,
            annotation.metas = "number2",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = genes),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
})




