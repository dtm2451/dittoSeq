# dittoSeq ![Logo](Vignette/dittoLogo_mini.png)

**A set of functions built to enable analysis and visualization of single-cell and bulk RNA-sequencing data by novice, experienced, and color blind coders**

dittoSeq includes universal plotting and helper functions for working with (sc)RNAseq data processed in these packages:

- Seurat (versions 2 & 3, single-cell RNAseq)
- SingleCellExperiment (single-cell RNAseq)
- DESeq2 (bulk RNAseq)
- edgeR (bulk RNAseq)
- Limma-Voom (bulk RNAseq)
- (Compatibility is planned for Monocle in a future version.)

All plotting functions spit out easy-to-read, color blind friendly, plots (ggplot2, plotly, or pheatmap) upon minimal coding input for your daily analysis needs, yet also allow sufficient manipulations to provide for out-of-the-box submission-quality figures!

dittoSeq also makes collection of underlying data easy, for submitting to journals, with `data.out = TRUE` inputs!

![Overview](Vignette/dittoSeq.gif)

Additionally, contains import functions for [Demuxlet](https://github.com/statgen/demuxlet) cell annotations as Mux-seq datasets often consist of side-by-side bulk and single-cell RNAseq.  (If you would like a pipeline for extraction of genotypes from bulk RNAseq to enable Demuxlet-calling of single-cell RNAseq, shoot me an email.)

## Installataion:

To install the version 0.3.0, run:

```
devtools::install_github("dtm2451/dittoSeq@Development-v0.3.0")
```

## News: version 0.3.0 is LIVE

**Includes lots of new features!**

  - Compatibility with bulk RNAseq data that was processed with edgeR & limma-voom.
  - `dittoScatterPlot()` which allows plotting of gene x gene / metadata x metadata / gene x metadata.  Great for examining raw droplet data QC, or potential marker gene RNA or CITE-seq expression.
  - `dittoDimPlot` now supports overlay of pseudotime-analysis trajectory paths.
  - Retrieval of underlying data as a dataframe or matrix.  Just add `data.out = TRUE` to the call.
  - Many more!

**Other updates since version 0.2:**

- Documentation overhaul for most functions with plenty of example code added.
- The package is now camelCase, **d**ittoSeq, to go along with typical conventions.
- Visualization names all start with `ditto...` and `multi_ditto...` instead of `DB...` and `multiDB...` Example: `dittoDimPlot` and `multi_dittoDimPlot`.
- Color Storage: The colors are now retrievable with a simple, empty, function call, `dittoColors()`.  For use outside of dittoSeq, simply use `dittoColors()`.  Example within dittoSeq: `dittoDimPlot("ident", seurat, color.panel = dittoColors() )`.

**More detailed vignettes are planned**, but in the meantime, if the in-R documentation (example: `?dittoScatterPlot`) and examples are not enough, please create an issue asking how to make whatever plot you would like!

## Color blindness friendliness:

The default colors of this package are meant to be color blind friendly.  To make it so, I used the suggested colors from this source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011](https://www.nature.com/articles/nmeth.1618) and adapted them slightly by appending darker and lighter versions to create a 24 color vector. All plotting functions use these colors, stored in `dittoColors()`, by default. Also included is a Simulate() function that allows you to see what your function might look like to a colorblind individual. For more info on that, see my [Colorblindness Compatibility Page](ColorblindCompatibility)

# Quick Start Guide:

(For rendered plots, download and open [Vignette/QuickStartRender.html](Vignette/QuickStartRender.html))

```{r}
# Install
devtools::install_github("dtm2451/dittoSeq@Development-v0.3.0")
# (Be sure to restart after a re-install!)
```

Load in your data, then go!:

```{r}
library(dittoSeq)
# library(Seurat)

# For working with scRNAseq data, works directly with Seurat and SingleCellExperiment objects
seurat <- Seurat::pbmc_small
dittoPlot("CD14", seurat, group.by = "ident")

sce <- Seurat::as.SingleCellExperiment(seurat)
dittoBarPlot("ident", sce, group.by = "RNA_snn_res.0.8")

# For working with bulk RNAseq data, first load your data into a format that dittoSeq quickly understands
# deseq2 <- importDESeq2()
# edger <- importEdgeR()
# limma.voom <- importEdgeR()
myRNA <- RNAseq_mock
dittoDimPlot("Gene1", myRNA, size = 3)
```

Quickly determine the metadata and gene options for plotting with helper functions:

```{r}
get.metas(seurat)
is.meta("nCount_RNA", seurat)

get.genes(myRNA)
is.gene("CD3E", myRNA)

# View them with these:
gene("CD3E", seurat, data.type = "raw")
meta("groups", seurat)
meta.levels("groups", seurat)
```

### There are many dittoSeq Plot Types

**Intuitive default adjustments generally allow creation of immediately useable plots.**

```{r}
# dittoPlot
dittoPlot("CD3E", seurat, group.by = "ident")
dittoPlot("CD3E", seurat, group.by = "ident",
    plots = c("boxplot", "jitter"))
dittoPlot("CD3E", seurat, group.by = "ident",
    plots = c("ridgeplot", "jitter"))

# dittoDimPlot
dittoDimPlot("ident", seurat, size = 3)
dittoDimPlot("CD3E", seurat, size = 3)

# dittoBarPlot
dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8")
dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    scale = "count")

# dittoHeatmap
dittoHeatmap(genes = get.genes(seurat)[1:20], seurat)
dittoHeatmap(genes = get.genes(seurat)[1:20], seurat,
    annotation.metas = c("groups", "ident"),
    scaled.to.max = TRUE,
    show.colnames = FALSE)
# Turning off cell clustering can be necessary for many cell scRNAseq
dittoHeatmap(genes = get.genes(seurat)[1:20], seurat,
    cluster_cols = FALSE)

# dittoScatterPlot
dittoScatterPlot(
    x.var = "CD3E", y.var = "nCount_RNA",
    color.var = "ident", shape.var = "RNA_snn_res.0.8",
    object = seurat,
    size = 3)

# Also:
    # multi_dittoDimPlot (multiple, in an array)
    # multi_dittoDimPlotVaryCells (multiple, in an array, but showing only certain
    #     cells in each plot)
    # multi_dittoPlot (multiple, in an array)
    # dittoPlot_VarsAcrossGroups (multiple genes or metadata as the jitterpoints (and
    #     other representations), summarized across groups by mean, median, ..., )
```

**Many adjustments can be made with simple additional inputs:**

Many adjustments to how data is reresented are within the examples above.  See documentation for more!  Also,

- All Titles are adjustable.
- Colors can be adjusted easily.
- Underlying data can be output.
- plotly hovering can be added.
- Many more! (Legends removal, label rotation, labels' and groupings' names, ...)

```{r}
dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    main = "Starters",
    sub = "By Type",
    xlab = NULL,
    ylab = "Generation 1",
    x.labels = c("Ash", "Misty"),
    legend.title = "Types",
    var.labels.rename = c("Fire", "Water", "Grass"),
    x.labels.rotate = FALSE)

dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    colors = c(3,1,2)) #Just changes the color order, probably most useful for dittoDimPlots
dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    color.panel = c("red", "orange", "purple"))

dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    data.out = TRUE)

dittoBarPlot("ident", seurat, group.by = "RNA_snn_res.0.8",
    do.hover = TRUE)
```
