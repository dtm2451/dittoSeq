# Tests for dittoPlotVarsAcrossGroups function
# library(dittoSeq); library(testthat); source("setup.R"); source("test-PlotVarsAcrossGroups.R")

pbmc@meta.data$number <- as.numeric(seq_along(colnames(pbmc)))
pbmc@meta.data$number2 <- as.numeric(seq_along(colnames(pbmc)))
genes <- get.genes(pbmc)[1:5]
grp <- "RNA_snn_res.1"
clr <- "orig.ident"
clr2 <- "RNA_snn_res.0.8"
cells.names <- colnames(pbmc)[1:40]
cells.logical <- c(rep(TRUE, 40), rep(FALSE,40))

test_that("dittoPlotVarsAcrossGroups can plot continuous metadata with all plot types", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter")),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("ridgeplot", "jitter")),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups can work for continuous metadata", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            c("number","number2"), pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter")),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups main legend can be removed", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            legend.show = FALSE),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups colors can be distinct from group.by", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            color.by = clr),
        "ggplot")
    expect_error(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            color.by = clr2),
        "Unable to interpret color.by input. All 'group.by' groupings must map to the same 'color.by' data.")
})

test_that("dittoPlotVarsAcrossGroups summary.fxn can be adjusted", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            summary.fxn = median),
        "ggplot")
    expect_error(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot", "jitter"),
            summary.fxn = function(x) x/2),
        NULL)
})

test_that("dittoPlotsVarsAcrossGroups can be subset to show only certain cells/samples with either cells.use method", {
    expect_s3_class(
        c1 <- dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = cells.names),
        "ggplot")
    expect_s3_class(
        c2 <- dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = cells.logical),
        "ggplot")
    expect_equal(c1,c2)
    # And if we remove an entire grouping...
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            cells.use = meta(grp,pbmc)!=0),
        "ggplot")
})

test_that("dittoPlots colors can be adjusted", {
    ### Manual check: These two should look the same.
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            color.panel = dittoColors()[5:1]),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            plots = c("vlnplot", "boxplot"),
            colors = 5:1),
        "ggplot")
})

test_that("dittoPlots titles and theme can be adjusted", {
    ### Manual check: All titles should be adjusted.
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            main = "Gotta catch", sub = "em all",
            xlab = "Pokemon", ylab = "Pokedex #s",
            legend.title = "groups"),
        "ggplot")
    ### Manual check: plot should be boxed
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            theme = theme_bw()),
        "ggplot")
})

test_that("dittoPlots y-axis can be adjusted, x for ridgeplots", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            min = -5, max = 100),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            y.breaks = seq(10,60,10)),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            min = -2, max = 1, plots = c("ridgeplot","jitter")),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            y.breaks = seq(-1,1.5,0.1), plots = c("ridgeplot","jitter")),
        "ggplot")
})

test_that("dittoPlots x-labels can be adjusted, (y) for ridgeplots", {
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            x.labels = 5:7),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            x.reorder = 3:1),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            x.labels.rotate = FALSE),
        "ggplot")
    ### Manual Check: L->R, green(5), blue(6), orange(7), with horizontal labels
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            x.labels = 5:7, x.reorder = 3:1, x.labels.rotate = FALSE),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups can have lines added", {
    # Manuel Check: Large blue dots that, in the yplot, look continuous accross groups.
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            add.line = 0),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            add.line = 0, line.linetype = "solid", line.color = "green"),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups jitter adjustments work", {
    # Manuel Check: Large blue dots that, in the yplot, look continuous accross groups.
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp, plots = "jitter",
            jitter.size = 10, jitter.color = "blue", jitter.width = 1),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups boxplot adjustments work", {
    # Manuel Check: Blue boxplots that touch eachother, with jitter visible behind in first plot, outliers shown in second
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp, plots = c("jitter", "boxplot"),
            boxplot.width = 1, boxplot.color = "blue", boxplot.fill = FALSE,
            boxplot.show.outliers = TRUE),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp, plots = c("jitter", "boxplot"),
            boxplot.width = 1, boxplot.color = "blue",
            boxplot.show.outliers = TRUE),
        "ggplot")
})

test_that("dittoPlotVarsAcrossGroups violin plot adjustments work", {
    # Manuel Check: Almost non-existent lines, with quite overlapping vlns.
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            vlnplot.lineweight = 0.1, vlnplot.width = 5),
        "ggplot")
    # The next three: first two look the same because equal numbers of dots, third should look different:
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            vlnplot.scaling = "count"),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            vlnplot.scaling = "area"),
        "ggplot")
    expect_s3_class(
        dittoPlotVarsAcrossGroups(
            genes, pbmc, group.by = grp,
            vlnplot.scaling = "width"),
        "ggplot")
})
