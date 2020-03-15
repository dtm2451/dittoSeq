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

test_that("Heatmap data type can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            slot = "counts"), # normally raw.data, but as.Seurat......
        "pheatmap")
})

# Set genes1:5 to have all zeros in logcounts in a new object
sce2 <- sce
assay(sce2,"logcounts")[1:5,] <- 0
test_that("Heatmap gives proper warnings when it should", {
    # Function throws a warning when any genes are not captured in the target cells.
    expect_warning(
        dittoHeatmap(
            genes,
            object = sce2,
            cells.use = meta("number", seurat)<5),
        "Gene(s) removed due to absence of expression within the 'cells.use' subset", fixed = TRUE)
    # Function throws an error when no genes provided are captured in the target cells.
    expect_error(
        dittoHeatmap(
            genes = genes[1:5],
            object = sce2,
            cells.use = meta("number", seurat)<10),
        "No target genes are expressed in the 'cells.use' subset", fixed = TRUE)
})

########################
##### Manual Check #####
########################

test_that("Heatmap title can be adjusted", {
    # Title chould be "Hello there!"
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            main = "Hello there!"),
        "pheatmap")
})

test_that("Heatmap sample renaming by metadata works", {
    ## Names should be numbers
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cell.names.meta = "number"),
        "pheatmap")
})

test_that("Heatmap can hide rownames/colnames", {
    # No colnames
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            show_colnames = FALSE),
        "pheatmap")
    # No rownames
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            show_rownames = FALSE),
        "pheatmap")
})

test_that("Heatmap highlight genes works", {
    # Only 1 label, gene1
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
    # Color bar should go from 0 to 1 and white to red
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            scaled.to.max = TRUE),
        "pheatmap")
})

test_that("Heatmap colors can be adjusted", {
    # yellow to black to red
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            heatmap.colors = colorRampPalette(c("yellow", "black", "red"))(50)),
        "pheatmap")
    # black to yellow, 0:1
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            scaled.to.max = TRUE,
            heatmap.colors.max.scaled = colorRampPalette(c("black", "yellow"))(25)),
        "pheatmap")
})

test_that("Heatmap annotations can be given & heatmaps can be ordered by metadata, expression, or user-input vector", {
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

test_that("Heatmap annotations can be given & ordering can be adjusted and follows defaults", {
    # Annotation bar = clusters
    # Cells should also be ordered by this (single-cell)
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "clusters"),
        "pheatmap")
    # Samples should not be ordered by this (bulk)
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = bulk,
            annotation.metas = "clusters"),
        "pheatmap")
    # Clusters even though an order.by would be given
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "clusters",
            cluster_cols = TRUE),
        "pheatmap")
    # Ordering, but distinct from the first annotation
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = c("clusters", "groups"),
            order.by = "groups"),
        "pheatmap")
})

test_that("Heatmap can be ordered when also subset to certain cells", {
    # Works for expression (ordered by gene1)
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            order.by = "gene1",
            cells.use = meta("number", seurat)<20),
        "pheatmap")
    # Works for metadata
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "groups",
            cells.use = colnames(seurat)[meta("number", seurat)<20]),
        "pheatmap")
    # Works with vectors provided
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number",
            order.by = seq_along(colnames(seurat)),
            cells.use = colnames(seurat)[meta("number", seurat)<20]),
        "pheatmap")
})

test_that("Heatmap can be subset to certain cells", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    # Few cells
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cells.use = meta("number", seurat)<10), # Logical method
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            cells.use = colnames(seurat)[meta("number", seurat)<10]), # names method
        "pheatmap")
    # Annotations still work (and ordering).
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "clusters",
            cells.use = colnames(seurat)[meta("number", seurat)<10]),
        "pheatmap")
})

test_that("Heatmap annotation colors can be adjusted", {
    # (via adjustment of the color pool)
    # green numeric (in order)
    # red, yellow, blue, purple for clusters
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = c("number","clusters"),
            annotation.colors = c("red", "yellow", "blue", "purple", "green3")),
        "pheatmap")
    # By annotation_colors (partial list!)
    color_list <- list(clusters = c('1' = "red",
                                    '2' = "yellow",
                                    '3' = "blue",
                                    '4' = "purple"))
    # Number color should change, but clusters should still be the same custom set!
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = c("number","clusters"),
            annotation_colors = color_list),
        "pheatmap")
})

test_that("Coloring works for discrete column and row annotations", {
    # column and row annotations, all discrete.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = c("clusters", "groups"),
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                genes,
                row.names = genes)),
        "pheatmap")
})

test_that("Coloring works for continuous column and row annotations", {
    # column and row annotations, all numeric.
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
    # column and row annotations, all numeric.
    # but column annotation bar flips, while legend stays the same.
    expect_s3_class(
        dittoHeatmap(
            genes,
            object = seurat,
            annotation.metas = "number2",
            order.by = "number",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = genes),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
})

test_that("scale and border_color pheatmap inputs function as expected", {
    expect_s3_class(
        dittoHeatmap(genes, object = seurat,
            scale = "column"),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(genes, object = seurat,
            scale = "none"),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(genes, object = seurat,
            scale = "column"),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(genes, object = seurat,
            border_color = "red"),
        "pheatmap")
})
