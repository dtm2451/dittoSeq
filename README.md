# dittoSeq

![Logo](Vignette/dittoLogo_mini.png)

**A set of functions built to enable analysis and visualization of single-cell and bulk RNA-sequencing data by novice, experienced, and color blind coders**

dittoSeq includes universal plotting and helper functions for working with (sc)RNAseq data processed in these packages:

- Seurat, (versions 2 & 3, single-cell RNAseq)
- SingleCellExperiment (single-cell RNAseq)
- DESeq2 (bulk RNAseq)
- (Compatibility is planned for Monocle, Limma-voom, and EdgeR)

All plotting functions spit out easy-to-read, color blind friendly, plots (ggplot2, plotly, or pheatmap) upon minimal coding input for your daily analysis needs, yet also allow sufficient manipulations to provide for out-of-the-box submission-quality figures!

dittoSeq also makes collection of underlying data easy, for submitting to journals, with `data.out = TRUE` inputs!

![Overview](Vignette/dittoSeq.gif)

Additionally, contains import functions for [Demuxlet](https://github.com/statgen/demuxlet) cell annotations as Mux-seq datasets often consist of side-by-side bulk and single-cell RNAseq.  (If you would like a pipeline for extraction of genotypes from bulk RNAseq to enable Demuxlet-calling of single-cell RNAseq, shoot me an email.)

## News: version 0.3.0 is on the way:

Includes lots of new features!

  - Bulk RNAseq compatibility for Limma-voom processed data [planning stages].
  - `dittoScatterPlot()` which allows plotting of gene x gene / metadata x metadata / gene x metadata.  Great for examining raw droplet data QC, or potential marker gene RNA or CITE-seq expression.
  - `dittoDimPlot` supports overlay of cluster-based trajectories paths in any dimensionality reduction space.
  - Retrieval of underlying data as a dataframe or matrix that is as simple as adding `data.out = TRUE`
  - Many more!

Other updates since version 0.2:

- Documentation overhaul for most functions with plenty of example code added.
- The package is now camelCase, **d**ittoSeq, to go along with typical conventions.
- Visualization names all start with `ditto...` and `multi_ditto...` instead of `DB...` and `multiDB...` Example: `dittoDimPlot` and `multi_dittoDimPlot`.
- Color Storage: The colors are now retrievable with a simple, empty, function call, `dittoColors()`.  For use outside of dittoSeq, simply use `dittoColors()`.  Example within dittoSeq: `dittoDimPlot("ident", seurat, color.panel = dittoColors() )`.

**Updated vignette with these functions is coming soon.** In the meantime, if the R documentation (example: `?dittoScatterPlot`) is not enough, please create an issue!

## To install:

To install the version 0.3.0, run:

```
devtools::install_github("dtm2451/dittoSeq@Development-v0.3.0")

```

## Color blindness friendliness:

The default colors of this package are meant to be color blind friendly.  To make it so, I used the suggested colors from this source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011](https://www.nature.com/articles/nmeth.1618) and adapted them slightly by appending darker and lighter versions to create a 24 color vector. All plotting functions use these colors, stored in `dittoColors()`, by default. Also included is a Simulate() function that allows you to see what your function might look like to a colorblind individual. For more info on that, see my [Colorblindness Compatibility Page](ColorblindCompatibility)
