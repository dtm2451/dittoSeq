# Tests for visualization functions
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Heatmap.R")

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
pbmc@meta.data$number2 <- as.numeric(rev(seq_along(colnames(pbmc))))

test_that("Coloring works for discrete column and row annotations", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_s3_class(
        dittoHeatmap(
            c("MS4A1","GNLY","CD3E","CD14","FCER1A",
              "FCGR3A","LYZ","PPBP","CD8A"),
            object = pbmc,
            col.annotation.metas = "ident",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                labels = c("MS4A1","GNLY","CD3E","CD14","FCER1A",
                    "FCGR3A","LYZ","PPBP","CD8A"),
                row.names = c("MS4A1","GNLY","CD3E","CD14","FCER1A",
                    "FCGR3A","LYZ","PPBP","CD8A"))),
        "pheatmap")
})

test_that("Coloring works for continuous column and row annotations", {
    # "Coloring works for discrete column annotations if BOTH PASS && LOOK DIFFERENT"
        # If: annotations are all continuous.
        # And...
    expect_s3_class(
        dittoHeatmap(
            c("MS4A1","GNLY","CD3E","CD14","FCER1A",
              "FCGR3A","LYZ","PPBP","CD8A"),
            object = pbmc,
            col.annotation.metas = "number",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = c("MS4A1","GNLY","CD3E","CD14","FCER1A",
                    "FCGR3A","LYZ","PPBP","CD8A")),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
        # And if column annotation bar flips, but legend stays the same.
    expect_s3_class(
        dittoHeatmap(
            c("MS4A1","GNLY","CD3E","CD14","FCER1A",
              "FCGR3A","LYZ","PPBP","CD8A"),
            object = pbmc,
            col.annotation.metas = "number2",
            scaled.to.max = TRUE,
            annotation_row = data.frame(
                lab = seq_len(9),
                row.names = c("MS4A1","GNLY","CD3E","CD14","FCER1A",
                    "FCGR3A","LYZ","PPBP","CD8A")),
            cluster_rows = FALSE,
            cluster_cols = FALSE),
        "pheatmap")
})

test_that("Heatmap highlight genes works", {
    # "Coloring works for discrete column and row annotations"
        # If: annotations are all discrete.
    expect_s3_class(
        dittoHeatmap(
            c("MS4A1","GNLY","CD3E","CD14","FCER1A",
              "FCGR3A","LYZ","PPBP","CD8A"),
            object = pbmc,
            highlight.genes = "CD3E"),
        "pheatmap")
})



