# dittoSeq Quick Start Guide:

![Logo](dittoLogo_mini.png)

For rendered plots, download and open [QuickStartRender.html](QuickStartRender.html)

*Note: dittoSeq also contains functions to aid in import of Demuxlet sample calls into Seurat single-cell RNAseq objects. For help with that, see the [Demuxlet Guide](../Demuxlet-Vignette/README.md).*

```{r}
# Install
devtools::install_github("dtm2451/dittoSeq")
# (Be sure to restart after a re-install!)
```

```{r}
# For older versions:
#   Old DB plotter version
# devtools::install_github("dtm2451/dittoSeq@v0.2")
#   Old DB plotter version plus some early ditto plotters
# devtools::install_github("dtm2451/dittoSeq@v0.2.20")
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
