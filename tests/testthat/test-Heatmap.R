# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Heatmap.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
seurat$number2 <- as.numeric(rev(seq_along(colnames(seurat))))
genes <- getGenes(seurat)[1:9]
metas <- c("score", "score2", "score3")

test_that("Heatmap can be plotted for Seurat or SCE", {
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = sce),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(
            metas = metas,
            object = sce),
        "pheatmap")
})

test_that("Heatmap data type can be adjusted", {
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            slot = "counts"), # normally raw.data, but as.Seurat......
        "pheatmap")
})

# Set genes1:5 to have all zeros in logcounts in a new object
sce2 <- sce
assay(sce2,"logcounts")[1:5,] <- 0

# Set metadata to zero as well
sce2$score <- 0
test_that("Heatmap gives warnings/errors when genes missing", {
    # Function throws a warning when any genes are not captured in the target cells.
    expect_warning(
        dittoHeatmap(
            genes = genes,
            object = sce2,
            cells.use = meta("number", seurat)<5),
        "Gene(s) or metadata removed due to absence of non-zero values within the 'cells.use' subset", fixed = TRUE)
    # Function throws an error when no genes provided are captured in the target cells.
    expect_error(
        dittoHeatmap(
            genes = genes[1:5],
            object = sce2,
            cells.use = meta("number", seurat)<10),
        "No target genes/metadata features have non-zero values in the 'cells.use' subset", fixed = TRUE)

    # And now for metadata.
    expect_warning(
        dittoHeatmap(
            genes = NULL,
            metas = metas,
            object = sce2,
            cells.use = meta("number", seurat)<5),
        "Gene(s) or metadata removed due to absence of non-zero values within the 'cells.use' subset", fixed = TRUE)
    # Function throws an error when no metas provided are captured in the target cells.
    expect_error(
        dittoHeatmap(
            genes = NULL,
            metas = metas[1],
            object = sce2,
            cells.use = meta("number", seurat)<10),
        "No target genes/metadata features have non-zero values in the 'cells.use' subset", fixed = TRUE)
})

test_that("Heatmap gives error when both genes and metas are not provided", {
    # Function throws an error when no genes or metas are provided.
    expect_error(
        dittoHeatmap(
            genes = NULL,
            metas = NULL,
            object = sce),
        "No 'genes' or 'metas' requested", fixed = TRUE)
})

test_that("Heatmap gives error when both highlight.genes and highlight.features are provided", {
    # Function throws an error when no genes or metas are provided.
    expect_error(
        dittoHeatmap(
            genes = genes,
            object = sce,
            highlight.features = "gene1",
            highlight.genes = "gene1"),
        "you can only specify one of 'highlight.genes' or 'highlight.features'", fixed = TRUE)
})

########################
##### Manual Check #####
########################

test_that("Heatmap title can be adjusted", {
    ### Title chould be "Hello there!"
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            main = "Hello there!"),
        "pheatmap")
})

test_that("Heatmap sample renaming by metadata works", {
    ### Names should be numbers
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            cell.names.meta = "number"),
        "pheatmap")
})

test_that("Heatmap can hide rownames/colnames", {
    ### No colnames
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            show_colnames = TRUE),
        "pheatmap")
    ### No rownames
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            show_rownames = FALSE),
        "pheatmap")
})

test_that("Heatmap highlight genes works", {
    ### Only 1 label, gene1
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            show_colnames = FALSE,
            show_rownames = FALSE,
            highlight.features = "gene1"),
        "pheatmap")
})

test_that("Heatmap can be scaled to max", {
    ### Color bar should go from 0 to 1 and white to red
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            scaled.to.max = TRUE),
        "pheatmap")
})

test_that("Heatmap colors can be adjusted", {
    ### yellow to black to red
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            heatmap.colors = colorRampPalette(c("yellow", "black", "red"))(50)),
        "pheatmap")
    ### black to yellow, 0:1
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            scaled.to.max = TRUE,
            heatmap.colors.max.scaled = colorRampPalette(c("black", "yellow"))(25)),
        "pheatmap")
})

test_that("Heatmap annotations can be given & heatmaps can be ordered by metadata, expression, or user-input vector", {
    # Works for expression
    ### ordered by gene1
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            order.by = "gene1"),
        "pheatmap")
    # Works for metadata
    ### ordered by groups
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            order.by = "groups",
            annot.by = "groups"),
        "pheatmap")
    # Works with vectors provided
    ### Ordered in REVERSE of number2
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "number2",
            order.by = seq_along(colnames(seurat))),
        "pheatmap")
})

test_that("Heatmap annotations can be given & ordering can be adjusted and follows defaults", {
    # Annotation bar = clusters
    ### Cells should also be ordered by this (single-cell)
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "clusters"),
        "pheatmap")
    # Samples should not be ordered by this (bulk)
    ### CLustered
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = bulk,
            annot.by = "clusters"),
        "pheatmap")
    # Clusters even though an order.by would be given
    ### Clustered!
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "clusters",
            cluster_cols = TRUE),
        "pheatmap")
    # Ordering, but distinct from the first annotation
    ### ordered by groups.
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = c("clusters", "groups"),
            order.by = "groups"),
        "pheatmap")
})

test_that("Heatmap can be ordered when also subset to certain cells", {
    # Works for expression
    ### Ordered by gene1
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            order.by = "gene1",
            cells.use = meta("number", seurat)<20),
        "pheatmap")
    # Works for metadata
    ### ordered by groups metadata
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "groups",
            cells.use = colnames(seurat)[meta("number", seurat)<20]),
        "pheatmap")
    # Works with vectors provided
    ### ordered in REVERSE of the number annotations
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "number2",
            order.by = seq_along(colnames(seurat)),
            cells.use = colnames(seurat)[meta("number", seurat)<20]),
        "pheatmap")
})

test_that("Heatmap can be subset to certain cells by any method", {
    # By logical
    ### Few cells
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            cells.use = meta("number", seurat)<10), # Logical method
        "pheatmap")
    # By names
    ### Same few cells
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            cells.use = colnames(seurat)[meta("number", seurat)<10]), # names method
        "pheatmap")
    # By indices
    ### Same few cells
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            cells.use = 1:9), # names method
        "pheatmap")
})

test_that("Heatmap annotation colors can be adjusted via annot.colors", {
    # (via adjustment of the color pool)
    ### red, yellow, blue, purple for clusters
    ### green numeric (in order)
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = c("number","clusters"),
            annot.colors = c("red", "yellow", "blue", "purple", "green3")),
        "pheatmap")
})

test_that("Heatmap annotation colors can be adjusted via annotation_colors", {
    color_list <- list(clusters = c('1' = "red",
                                    '2' = "yellow",
                                    '3' = "blue",
                                    '4' = "purple"))
    # Number color should change, but clusters should still be the same custom set as above.
    ### red, yellow, blue, purple for clusters
    ### dittoBlue for number
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = c("number","clusters"),
            annotation_colors = color_list),
        "pheatmap")
    ### When annotation_colors provides all colors for annotation_col.
    color_list <- list(clusters = c('1' = "red",
                                    '2' = "yellow",
                                    '3' = "blue",
                                    '4' = "purple"))
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = c("clusters"),
            annotation_colors = color_list),
        "pheatmap")
})

test_that("Coloring works for discrete column and row annotations", {
    ### column and row annotations, all discrete & all dittoColors
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = c("clusters", "groups"),
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                genes,
                row.names = genes)),
        "pheatmap")
})

test_that("Coloring works for continuous column and row annotations", {
    ### column and row annotations, all numeric.
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "number",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = genes),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
    ### column and row annotations, all numeric.
    ### but column annotation bar flips, while legend stays the same.
    expect_s3_class(
        dittoHeatmap(
            genes = genes,
            object = seurat,
            annot.by = "number2",
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
        dittoHeatmap(genes = genes, object = seurat,
            scale = "none"),
        "pheatmap")
    expect_s3_class(
        dittoHeatmap(genes = genes, object = seurat,
            border_color = "red"),
        "pheatmap")
})

test_that("dittoHeatmap swap.rownames works", {
    swap_genes <- paste(genes, "symb", sep = "_")
    
    expect_s3_class(
        dittoHeatmap(genes = swap_genes, object = sce, swap.rownames = "symbol"),
        "pheatmap")
    expect_equivalent(
        rownames(
            dittoHeatmap(
                genes = swap_genes, object = sce, swap.rownames = "symbol",
                data.out = TRUE
            )$mat),
        swap_genes)
})
