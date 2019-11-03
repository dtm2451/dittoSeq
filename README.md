# dittoSeq

![Logo](Vignette/dittoLogo_mini.png)

**A set of functions built to enable analysis and visualization of single-cell and bulk RNA-sequencing data by novice, experienced, and color blind coders**

dittoSeq includes universal plotting and helper functions for working with (sc)RNAseq data processed in these packages:

- Seurat, (versions 2 & 3, single-cell RNAseq)
- SingleCellExperiment (single-cell RNAseq)
- DESeq2 (bulk RNAseq)
- edgeR (bulk RNAseq)
- Limma-Voom (bulk RNAseq)
- (Compatibility is planned for Monocle in a future version.)

All plotting functions spit out easy-to-read, color blind friendly, plots (ggplot2, plotly, or pheatmap) upon minimal coding input for your daily analysis needs, yet also allow sufficient manipulations to provide for out-of-the-box submission-quality figures!

dittoSeq also makes collection of underlying data easy, for submitting to journals, with `data.out = TRUE` inputs!

![Overview](Vignette/dittoSeq.gif)

Additionally, contains import functions for [Demuxlet](https://github.com/statgen/demuxlet) cell annotations as Mux-seq datasets often consist of side-by-side bulk and single-cell RNAseq.  (If you would like a pipeline for extraction of genotypes from bulk RNAseq to enable Demuxlet-calling of single-cell RNAseq, shoot me an email.)

## News: version 0.3.0 is on the way:

Includes lots of new features!

  - Compatibility with bulk RNAseq data that was processed with edgeR & limma-voom.
  - `dittoScatterPlot()` which allows plotting of gene x gene / metadata x metadata / gene x metadata.  Great for examining raw droplet data QC, or potential marker gene RNA or CITE-seq expression.
  - `dittoDimPlot` now supports overlay of pseudotime-analysis trajectory paths.
  - Retrieval of underlying data as a dataframe or matrix.  Just add `data.out = TRUE` to the call.
  - Many more!

Other updates since version 0.2:

- Documentation overhaul for most functions with plenty of example code added.
- The package is now camelCase, **d**ittoSeq, to go along with typical conventions.
- Visualization names all start with `ditto...` and `multi_ditto...` instead of `DB...` and `multiDB...` Example: `dittoDimPlot` and `multi_dittoDimPlot`.
- Color Storage: The colors are now retrievable with a simple, empty, function call, `dittoColors()`.  For use outside of dittoSeq, simply use `dittoColors()`.  Example within dittoSeq: `dittoDimPlot("ident", seurat, color.panel = dittoColors() )`.

**Updated vignette with these functions is coming soon.** In the meantime, if the R documentation (example: `?dittoScatterPlot`) is not enough, please create an issue!  For working with bulk data, start with `?importDESeq2` if the data was analyzed with DESeq2 or `?importEdgeR` if the data was analyzed with edgeR or Limma-Voom.

A quick additional note on use of dittoSeq for bulkRNAseq data: dittoSeq is built to be a companion package for aiding analysis of RNAseq datasets of diverse structures. It works quite well for this purpose, but is not meant to be used as the main tool for data-processing and normalization. Inclusion of automatic normalization of bulk RNAseq data in import functions is provided for ease-of-use, but it is the responsibility of the user to provide properly normalized data (to the `normalized.data` inputs of import functions) when the default settings for rlog (DESeq2) or cpm (edgeR/Limma-Voom) calculation are improper.

## To install:

To install the version 0.3.0, run:

```
devtools::install_github("dtm2451/dittoSeq@Development-v0.3.0")
```

## Color blindness friendliness:

The default colors of this package are meant to be color blind friendly.  To make it so, I used the suggested colors from this source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011](https://www.nature.com/articles/nmeth.1618) and adapted them slightly by appending darker and lighter versions to create a 24 color vector. All plotting functions use these colors, stored in `dittoColors()`, by default. Also included is a Simulate() function that allows you to see what your function might look like to a colorblind individual. For more info on that, see my [Colorblindness Compatibility Page](ColorblindCompatibility)
