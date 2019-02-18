# DittoSeq

**A set of functions built to enable analysis and visualization of single cell and bulk RNA-sequencing data by novice and/or color blind coders**

![example1](Vignette/DBDimPlot2.png)
![example2](Vignette/DBBarPlot1.png)
![example3](Vignette/DBPlot2.png)
![example4](Vignette/multiDBDimPlot1.png)

*For a tutorial on how to use these functions, click [here](Vignette)*

Package includes various helper and plotting functions for working with RNAseq data analyzed in other packages; Seurat for single-cell RNAseq data, DESeq for bulk RNAseq data.

All plotting functions spit out easy-to-read, color blind friendly, ggplot plots upon minimal coding input for your daily analysis needs, and they also allow sufficient manipulations to provide for out-of-the-box submission-quality figures.

In addition to significantly extending visualization functionality of the widely used Seurat package for single cell RNAseq data, the package allows generation of similar figures from bulk RNA sequencing data. Thus, it enables analysis of single cell and bulk data side-by-side.

NOTE: I use this package daily, and am constantly coming up with new ideas for tweaks and additional utility myself.  To report errors, give feedback, or suggest new features, you can do so either through [github](https://github.com/dtm2451/DittoSeq), or by email at <daniel.bunis@ucsf.edu>.

## Color blindness friendliness:

The default colors of this package are meant to be color blind friendly.  To make it so, I used the suggested colors from this source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011.](https://www.nature.com/articles/nmeth.1618)  Darker and lighter versions of these same colors were then appended to make it a 24 color vector.  All plotting functions use these colors, stored in `MYcolors`, by default.  Also included is a `Simulate()` function that allows you to see what your function might look like to a colorblind individual.  For more info on that, see my [Colorblindness Compatibility Page](ColorblindCompatibility)

## To use:

Simply run this code.  It may take a while due to installation of Seurat and DESeq.  I may remove the need for Seurat in a future version:

```
devtools::install_github("dtm2451/DittoSeq")
```

For an explanation on how to use these functions, see [my Vignette](Vignette)

## Plotting Functions

**`DBDimPlot()`** = handles all needs for Seurat TSNEPlot / PCAPlot / DimPlot functions.  Improves on the Seurat functions' capabilities to present continuous (including negative) numerical data, or descrete data (clustering, samples, batches, condition, etc.) in various ways.

**`DBPlot()`** = handles needs of Seurat's VlnPlot function. Allows generation of jitter/dot-plot, boxplot, and/or violin-plot representation of numerical data, with order of what's on top easily settable. Data can be expression of particular genes or any numerical metadata like percent.mito, nUMI, and nGene.  Colors and grouping of cells is tunable through discrete inputs.

**`DBBarPlot()`** = No analogous function currently in Seurat, which is a bit crazy imho. Most common use: Plotting the cluster breakdown of all cells of each sample. Essentially, it is similar to DBPlot, but for discrete variables. Handles plotting of discrete data on a per-sample or per-condition grouping.

**multi-plotters** = Plot multiple DBDimPlots or DBPlots in an array.  Can handle most inputs that would be given to the individual functions.  Names are **`multiDBDimPlot()`**, **`multiDBPlot()`**, and **`multiDBDimPlot_vary_cells()`**.

## Helper functions

These make manipulating Seurat data, and using my plotting functons, easier.

**`get.metas()`** and **`get.genes()`**: Returns the list of meta.data slots or the list of genes included in the dataset.  Works exactly like typing `names(object@meta.data)` or `rownames(object@raw.data)`, only easier.

**`is.meta()`** and **`is.gene()`**: Returns TRUE or FALSE for whether a "meta.data" or "gene" input is part of the dataset.

**`meta()`** and **`gene()`**: Returns the values of a meta.data for every cell or the normalized expression data (`@data` slot) for all cells.

**`meta.levels()`**: Returns the range of values of metadata. Like running `levels(as.factor(object@meta.data$meta))`. Alternatively, can reurn the counts of each value of the meta.data if the optional input `table.out` is set to `TRUE`.

**`extDim()`**: extracts the loadings of each cell for a given dimensional reduction space.

## Color adjustment functions

**`Darken()`**: Darkens a color or color.panel by a given amount. (note: use these on a color.panel, not on a generated plot)

**`Lighten()`**: Lightens a color or color.panel by a given amount. (note: use these on a color.panel, not on a generated plot)

**`Simulate()`**: Generates any of the plot-types included in this package with colors adjusted to simulate any of the major forms of colorblindness.
