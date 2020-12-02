# Tests for dittoPlot function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-Plot.R")

seurat$number <- as.numeric(seq_along(colnames(seurat)))
grp <- "clusters"
clr <- "age"
clr2 <- "groups"
gene1 <- "gene1"
cells.names <- colnames(seurat)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,ncells-40))

test_that("dittoPlot can plot continuous metadata with all plot types", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter")),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("ridgeplot", "jitter")),
        "ggplot")
})

test_that("dittoPlot can work for SCE", {
    expect_s3_class(
        dittoPlot(
            gene1, object=sce, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter")),
        "ggplot")
})

test_that("dittoPlot can plot gene expression data with all plot types", {
    expect_s3_class(
        dittoPlot(
            gene1, object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter")),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            gene1, object=seurat, group.by = grp,
            plots = c("ridgeplot", "jitter")),
        "ggplot")
})

test_that("dittoPlots can be subset to show only certain cells/samples with any cells.use method", {
    expect_s3_class(
        c1 <- dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = cells.logical),
        "ggplot")
    expect_s3_class(
        c3 <- dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = 1:40),
        "ggplot")
    expect_equal(c1,c2)
    expect_equal(c1,c3)
    # And if we remove an entire grouping...
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = meta(grp,seurat)!=0),
        "ggplot")
})

test_that("dittoPlot main legend can be removed", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            legend.show = FALSE),
        "ggplot")
})

test_that("dittoPlot colors can be distinct from group.by", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            color.by = clr),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            color.by = clr2),
        "ggplot")
})

test_that("dittoPlot shapes can be a metadata and distinct from group.by", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.by = grp),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.by = clr2),
        "ggplot")
})

test_that("dittoPlot shapes can be adjusted in many ways", {
    # Shapes should be triangles instead of dots
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.panel = 17),
        "ggplot")
    # Shapes should be dot and triangle instead of dot and square
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.by = clr2, shape.panel = 16:19),
        "ggplot")
    # Shapes should be enlarged in the legend
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.by = clr2, jitter.shape.legend.size = 5),
        "ggplot")
    # Shapes legend should be removed
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            shape.by = clr2, jitter.shape.legend.show = FALSE),
        "ggplot")
})

test_that("dittoPlots colors can be adjusted", {
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            color.panel = dittoColors()[5:1]),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            colors = 5:1),
        "ggplot")
})

test_that("dittoPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: plot should be boxed
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoPlots y-axis can be adjusted, x for ridgeplots", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            min = -5, max = 100),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            y.breaks = seq(10,60,10)),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            min = -50, max = 100),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            y.breaks = seq(10,60,10)),
        "ggplot")
})

test_that("dittoPlots x-labels can be adjusted, (y) for ridgeplots", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            x.labels = 5:8),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            x.reorder = 4:1),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            x.labels.rotate = FALSE),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            x.labels = 5:8),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            x.reorder = 4:1),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: L->R, green(5), blue(6), orange(7), with horizontal labels
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            x.labels = 5:8, x.reorder = 4:1, x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: B -> T, green(5), blue(6), orange(7), with rotated labels
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            x.labels = 5:8, x.reorder = 4:1, x.labels.rotate = TRUE),
        "ggplot")
})

test_that("dittoPlot can have lines added", {
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            add.line = 20),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            add.line = 20, line.linetype = "solid", line.color = "green"),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp,
            add.line = 20, line.linetype = "solid", line.color = "green"),
        "ggplot")
})

test_that("dittoPlot jitter adjustments work", {
    # Manuel Check: Large blue dots that, in the yplot, look continuous accross groups.
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp, plots = "jitter",
            jitter.size = 10, jitter.color = "blue", jitter.width = 1),
        "ggplot")
    expect_s3_class(
        dittoRidgePlot(
            "number", object=seurat, group.by = grp, plots = c("jitter","ridgeplot"),
            jitter.size = 10, jitter.color = "blue", jitter.width = 1),
        "ggplot")
})

test_that("dittoPlot boxplot adjustments work", {
    # Manuel Check: Blue boxplots that touch eachother, with jitter visible behind.
    # Not actually checked here manually: whether outliers are shown cuz there are none.
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp, plots = c("jitter", "boxplot"),
            boxplot.width = 1, boxplot.color = "blue", boxplot.fill = FALSE,
            boxplot.show.outliers = TRUE),
        "ggplot")
    # Manual Check: boxplots that overlap
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp, plots = c("jitter","boxplot"),
            color.by = clr, boxplot.width = 0.3, boxplot.position.dodge = 0.2),
        "ggplot")
})

test_that("dittoPlot violin plot adjustments work", {
    # Manuel Check: Almost non-existent lines, with quite overlapping vlns.
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            vlnplot.lineweight = 0.1, vlnplot.width = 5),
        "ggplot")
    # The next three should look different from eachother:
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            vlnplot.scaling = "count"),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            vlnplot.scaling = "area"),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            "number", object=seurat, group.by = grp,
            vlnplot.scaling = "width"),
        "ggplot")
})

test_that("dittoPlot can add extra vars to dataframe", {
    df1 <- dittoPlot(
            gene1, grp, grp, object = seurat,
            data.out = TRUE)[[2]]
    expect_s3_class(
        df2 <- dittoPlot(
            gene1, grp, grp, object = seurat,
            extra.vars = c(clr, clr2), data.out = TRUE)[[2]],
        "data.frame")
    expect_equal(ncol(df1), 3)
    expect_equal(ncol(df2), 5)
})

test_that("dittoPlot can be facted with split.by (1 or 2 vars)", {
    # MANUAL CHECK: FACETING
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            split.by = grp),
        "ggplot")
    # horizontal
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            split.by = grp,
            split.nrow = 1),
        "ggplot")
    # vertical
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            split.by = grp,
            split.ncol = 1),
        "ggplot")
    # Grid with rows=age, cols=groups
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            split.by = c(grp,clr)),
        "ggplot")
    # Works with cells.use (should have grey cells)
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            split.by = c(grp,clr),
            cells.use = cells.logical),
        "ggplot")
})

test_that("dittoPlot with and without jitter rasterization produces identical plots", {
    # MANUAL CHECK: Should be indentical
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            do.raster = TRUE),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat),
        "ggplot")
    # Jitter in front.
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            do.raster = TRUE, 
            plots = c("vlnplot", "jitter")),
        "ggplot")
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat, 
            plots = c("vlnplot", "jitter")),
        "ggplot")
    # Should be lower resolution
    expect_s3_class(
        dittoPlot(
            gene1, grp, grp, object = seurat,
            do.raster = TRUE,
            raster.dpi = 10,
            plots = c("vlnplot", "jitter")),
        "ggplot")
})

test_that("dittoPlot swap.rownames works", {
    expect_s3_class(
        dittoPlot(sce, "gene1_symb", grp, swap.rownames = "symbol"),
        "ggplot")
})
