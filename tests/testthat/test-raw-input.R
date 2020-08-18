# setwd("dittoSeq/tests/testthat")
# library(dittoSeq); library(testthat); source("setup.R"); source("test-raw-input.R")

# For data targeting
sce$number <- seq_len(ncol(sce))
disc1 <- "groups"
disc2 <- "clusters"
disc3 <- "age"
cont1 <- "number"
cont2 <- "score"
cont3 <- "score2"
cont4 <- "score3"
raw <- as.vector(gene("gene1", sce))

# For cell/sample subsetting
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))

# For checking different data structures
colnames(sce) <- paste0("cell", seq_len(ncol(sce)))

# Extract raw
count_mtx <- counts(sce)
meta_DF <- colData(sce)
meta_df <- data.frame(meta_DF)
pc_mtx <- data.frame(reducedDim(sce))

#####################
### Base Function ###
#####################

test_that("raw object creation makes sce with all parts", {
    expect_s4_class(
        x <- dittoSeq:::.make_sce_if_raw(
            object = count_mtx,
            metadata = meta_DF,
            reduction.matrix = pc_mtx),
        "SingleCellExperiment")
    
    expect_equal(counts(x), count_mtx)
    expect_equal(getMetas(x, FALSE), meta_DF)
    expect_equal(dittoSeq:::.extract_Reduced_Dim("Dim", 1, x)$embeddings, pc_mtx[,1])
})

test_that("raw object creation makes sce when given all dataframes", {
    expect_s4_class(
        x <- dittoSeq:::.make_sce_if_raw(
            object = data.frame(count_mtx),
            metadata = meta_df,
            reduction.matrix = data.frame(pc_mtx)),
        "SingleCellExperiment")
    
    expect_equal(as.matrix(counts(x)), as.matrix(count_mtx))
    expect_equal(getMetas(x, FALSE), meta_DF)
    expect_equal(dittoSeq:::.extract_Reduced_Dim("Dim", 1, x)$embeddings, pc_mtx[,1])
})

test_that("raw object creation makes sce when given all DataFrames", {
    expect_s4_class(
        x <- dittoSeq:::.make_sce_if_raw(
            object = DataFrame(count_mtx),
            metadata = meta_DF,
            reduction.matrix = DataFrame(pc_mtx)),
        "SingleCellExperiment")
    
    expect_equal(as.matrix(counts(x)), as.matrix(count_mtx))
    expect_equal(getMetas(x, FALSE), meta_DF)
    expect_equal(dittoSeq:::.extract_Reduced_Dim("Dim", 1, x)$embeddings, pc_mtx[,1])
})

test_that("raw object creation makes sce when numeric data given as DelayedArrays", {
    expect_s4_class(
        x <- dittoSeq:::.make_sce_if_raw(
            object = DelayedArray(count_mtx),
            metadata = meta_DF,
            reduction.matrix = DelayedArray(pc_mtx)),
        "SingleCellExperiment")
    
    expect_equal(as.matrix(counts(x)), as.matrix(count_mtx))
    expect_equal(getMetas(x, FALSE), meta_DF)
    expect_equal(dittoSeq:::.extract_Reduced_Dim("Dim", 1, x)$embeddings, pc_mtx[,1])
})

##########################
### Plot Compatibility ###
##########################

test_that("raw input works for dittoScatterPlot", {
    expect_s3_class(
        dittoScatterPlot(
            "gene1", "gene2",
            object = count_mtx),
        "ggplot")
    
    expect_s3_class(
        dittoScatterPlot(
            cont1, cont2, disc1,
            cells.use = cells.names,
            adjustment.x = "z-score",
            # rename.color.groups = 1:5, order = "decreasing",
            do.label = TRUE, labels.highlight = FALSE, labels.repel = FALSE,
            add.trajectory.lineages = list(
                c("B","A","C"),
                c("C","A")),
            trajectory.cluster.meta = disc1,
            object = count_mtx,
            metadata = meta_df),
        "ggplot")
})

test_that("raw input works for dittoDimPlot", {
    # Single
    expect_s3_class(
        dittoDimPlot(
            "gene1",
            object = count_mtx,
            reduction.use = pc_mtx),
        "ggplot")
    
    # Multi
    expect_s3_class(
        multi_dittoDimPlot(
            c("gene1", "gene2", disc1),
            object = count_mtx,
            metadata = meta_df,
            reduction.use = pc_mtx),
        "gtable")
    
    # Single - Complicated
    expect_s3_class(
        dittoDimPlot(
            disc1,
            cells.use = cells.names,
            rename.var.groups = 1:5, order = "decreasing",
            size = 8,
            do.label = TRUE, labels.highlight = FALSE, labels.repel = FALSE,
            do.contour = TRUE,
            add.trajectory.curves = list(
                data.frame(
                    c(-1,0,-2),
                    c(-2,-1,0)),
                data.frame(
                    c(-1:1),
                    c(2,0,1)
                )),
            add.trajectory.lineages = list(
                c("B","A","C"),
                c("C","A")),
            trajectory.cluster.meta = disc1,
            object = count_mtx,
            reduction.use = pc_mtx,
            metadata = meta_df,
            split.by = "age"),
        "ggplot")
    expect_s3_class(
        dittoDimPlot(
            disc1,
            cells.use = cells.names,
            size = 8,
            do.ellipse = TRUE,
            do.letter = TRUE,
            object = count_mtx,
            reduction.use = pc_mtx,
            metadata = meta_df),
        "ggplot")
})

test_that("raw input works for dittoScatterHex", {
    expect_s3_class(
        dittoScatterHex(
            "gene1", "gene2",
            object = count_mtx),
        "ggplot")
})

test_that("raw input works for dittoDimHex", {
    expect_s3_class(
        dittoDimHex(
            object = count_mtx,
            reduction.use = pc_mtx),
        "ggplot")
})

test_that("raw input works for dittoPlot", {
    # Single
    expect_s3_class(
        dittoPlot(
            "gene1", disc1,
            object = count_mtx,
            metadata = meta_df),
        "ggplot")
    
    # Multi
    expect_s3_class(
        multi_dittoPlot(
            c("gene1", "gene2"), disc1,
            object = count_mtx,
            metadata = meta_df),
        "gtable")
    
    # Single - Complicated
    expect_s3_class(
        dittoPlot(
            "gene1", disc1,
            color.by = disc2,
            shape.by = disc3,
            split.by = disc2,
            cells.use = 1:70,
            adjustment = "relative.to.max",
            object = count_mtx,
            metadata = meta_df),
        "ggplot")
})

test_that("raw input works for dittoBarPlot", {
    expect_s3_class(
        dittoBarPlot(
            disc1, disc2,
            object = count_mtx,
            cells.use = 1:10,
            scale = "count",
            metadata = meta_df),
        "ggplot")
})

test_that("raw input works for dittoHeatmap", {
    expect_s3_class(
        dittoHeatmap(
            object = count_mtx),
        "pheatmap")
    
    # More complicated
    expect_s3_class(
        dittoHeatmap(
            cells.use = 1:10,
            annot.by = c(disc1, cont1),
            cell.names.meta = "number",
            object = count_mtx,
            metadata = meta_df),
        "pheatmap")
})

test_that("raw input works for dittoPlotVarsAcrossGroups", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            object = count_mtx,
            vars = c("gene1", "gene2", "gene3", "gene4", "gene5"),
            group.by = disc1,
            cells.use = 1:10, 
            metadata = meta_df),
        "ggplot")
})

test_that("raw input works for multi_dittoDimPlotVaryCells", {
    expect_s3_class(
        multi_dittoDimPlotVaryCells(
            cont1,
            vary.cells.meta = disc1,
            object = count_mtx,
            reduction.use = pc_mtx,
            metadata = meta_df),
        "gtable")
    expect_s3_class(
        multi_dittoDimPlotVaryCells(
            disc1,
            vary.cells.meta = disc1,
            object = count_mtx,
            reduction.use = pc_mtx,
            metadata = meta_df),
        "gtable")
})

######################
### Error Messages ###
######################

test_that("meta() will give proper error when no metadata is provided", {
    expect_error(dittoPlot(
        "gene1", object = count_mtx,
        group.by = "a"),
        "no metadata", fixed = TRUE)
})

test_that("meta() messages when original colnames(object) not used", {
    count_mtx_noNames <- count_mtx
    colnames(count_mtx_noNames) <- NULL
    expect_message(dittoPlot(
        "gene1", object = count_mtx_noNames, metadata = meta_df,
        group.by = "groups"),
        "NULL or invalid, using rownames('metadata').", fixed = TRUE)
    expect_message(dittoHeatmap(object = count_mtx_noNames),
        "NULL or invalid, using col1, col2, col3, etc.", fixed = TRUE)
})
