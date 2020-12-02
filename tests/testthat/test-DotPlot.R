# Tests for dittoDotPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-DotPlot.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
seurat$number2 <- rev(seurat$number)
disc <- "clusters"
disc2 <- "timepoint"
disc3 <- "age"
genes <- getGenes(seurat)[1:5]
metas <- c("score", "score2", "score3")
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))

test_that("dittoDotPlot can plot gene and meta data with Seurat or SCE", {
    expect_s3_class(
        print(dittoDotPlot(seurat, group.by = disc,
            genes)),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(sce, group.by = disc,
            genes),
        "ggplot")
    
    expect_s3_class(
        dittoDotPlot(seurat, group.by = disc,
            metas),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(sce, group.by = disc,
            metas),
        "ggplot")
    
    expect_s3_class(
        dittoDotPlot(seurat, group.by = disc,
            c("score", "gene1")),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(sce, group.by = disc,
            c("score", "gene1")),
        "ggplot")
})

test_that("dittoDotPlot errors for single vars or non-numeric vars", {
    expect_error(
        dittoDotPlot(seurat, group.by = disc,
            c("score")),
        "'vars' must be a vector of at least two", fixed = TRUE)
    
    expect_error(
        dittoDotPlot(seurat, group.by = disc,
            c("gene1")),
        "'vars' must be a vector of at least two", fixed = TRUE)
    expect_error(
        dittoDotPlot(seurat, group.by = disc,
            c("gene1", "gene2", disc)),
        "'vars' must be numeric", fixed = TRUE)
})

test_that("dittoDotPlot works with any gene adjustments", {
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            adjustment = "relative.to.max"),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            adjustment = "z-score"),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            adjustment = NULL),
        "ggplot")
})

test_that("dittoDotPlot data.out works", {
    # scaling on by default
    expect_type(
        d <- dittoDotPlot(seurat, genes, disc,
            data.out = TRUE),
        "list")
    
    expect_equal(names(d), c("p", "data"))
    expect_s3_class(d$p, "ggplot")
    expect_s3_class(d$data, "data.frame")
})

test_that("dittoDotPlot scaling acts as expected", {
    # scaling on by default
    expect_type(
        d <- dittoDotPlot(seurat, genes, disc,
            data.out = TRUE),
        "list")
    # Controlled by 'scale'
    expect_type(
        d2 <- dittoDotPlot(seurat, genes, disc,
            scale = FALSE,
            data.out = TRUE),
        "list")
    
    expect_true( identical(d2$data$color, d$data$pre.scale))
    expect_false(identical(d2$data$color, d$data$color))
})

test_that("dittoDotPlot legend can be removed", {
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            legend.show = FALSE),
        "ggplot")
})

test_that("dittoDotPlot summary.fxns can be adjusted", {
    # Manual check: scale of color in 1, should have same range as size in 2.
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            summary.fxn.color = median),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            summary.fxn.size = median),
        "ggplot")
    
    expect_error(
        dittoDotPlot(seurat, genes, disc,
            summary.fxn.color = function(x) x/2),
        "result is length", fixed = TRUE)
})

test_that("dittoDotPlot colors, sizes, ranges, legends are adjustable", {
    # Manual check: color = black to light grey, -3 to 3.
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            min.color = "black", max.color = "grey90",
            min = -3, max = 3),
        "ggplot")
    
    # Manual check: size = -0.25 to 1, dots LARGE
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            min.percent = -0.25, max.percent = 1,
            size = 10),
        "ggplot")
    
    # Manual check: color = 3 breaks, low, 0, high
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            min = -3, max = 3,
            legend.color.breaks = seq(-3,3,3),
            legend.color.breaks.labels = c("low", 0, "high")),
        "ggplot")
})

test_that("dittoDotPlot can be subset to show only certain cells/samples with any cells.use method", {
    expect_s3_class(
        c1 <- dittoDotPlot(seurat, genes, disc,
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoDotPlot(seurat, genes, disc,
            cells.use = cells.logical),
        "ggplot")
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            cells.use = 1:40),
        "ggplot")
    expect_equal(c1,c2)
    # And if we remove an entire grouping...
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            cells.use = meta(disc,seurat)!=0),
        "ggplot")
})

test_that("dittoDotPlot titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.size.title = "Pika",
            legend.color.title = "chu"),
        "ggplot")
    ### Manual check: plot should be boxed with grid lines
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoDotPlot y-labels can be adjusted", {
    # Manual check: Labels = 3,4,5,6
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            y.labels = 3:6),
        "ggplot")
    # Manual check: 4 at bottom
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            y.reorder = 4:1),
        "ggplot")
    # Manual check: x-labels horizontal
    expect_s3_class(
        dittoDotPlot(seurat, genes, disc,
            x.labels.rotate = FALSE),
        "ggplot")
})

test_that("dittoDotPlot assay, slot, adjustment work", {
    expect_type(
        d_raw <- dittoDotPlot(seurat, genes, disc, data.out = TRUE, scale = FALSE,
            slot = "counts"),
        "list")
    expect_type(
        d_log <- dittoDotPlot(seurat, genes, disc, data.out = TRUE, scale = FALSE,
            slot = "data"),
        "list")
    expect_true(all(
        d_raw$data$color >= d_log$data$color))
    
    expect_type(
        d_raw <- dittoDotPlot(sce, genes, disc, data.out = TRUE, scale = FALSE,
            assay = "counts"),
        "list")
    expect_type(
        d_log <- dittoDotPlot(sce, genes, disc, data.out = TRUE, scale = FALSE,
            slot = "logcounts"),
        "list")
    expect_true(all(
        d_raw$data$color >= d_log$data$color))
})

test_that("dittoDotPlot swap.rownames works", {
    
    swap_genes <- paste(genes, "symb", sep = "_")
    
    no_swap <- dittoDotPlot(sce, genes, disc, data.out = TRUE)
    swap <- dittoDotPlot(sce, swap_genes, disc, data.out = TRUE,
                         swap.rownames = "symbol")
    
    expect_equivalent(no_swap$data$color, swap$data$color)
    expect_equivalent(swap$data$var,
                      factor(paste(no_swap$data$var, "symb", sep = "_")))
    
    expect_s3_class(
        no_swap$p,
        "ggplot")
    expect_s3_class(
        swap$p,
        "ggplot")
})
