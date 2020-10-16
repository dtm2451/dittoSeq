\name{dittoSeqNews}
\title{dittoSeq News}
\encoding{UTF-8}

\section{Version 1.1.12}{\itemize{
\item `dittoDotPlot()` `vars` ordering is now set based on the order given (adds user control compared the the previous alphabetical ordering.)
\item Bug fix: `do.hover` in `dittoDotPlot()` now shows actual hover info rather than "hover.string".
\item Bug fix: Corrected the alignment mismatch for `dittoPlot()` box and violin plots when `group.by` != `color.by`
\item Vignette & documentation updates
}}

\section{Version 1.1.11}{\itemize{
\item Added compatibility for SummarizedExperiments natively, rather than requiring their conversion to a SingleCellExperiment.
\item Added rasterization & ComplexHeatmap compatibility to `dittoHeatmap()` via utilization of the ComplexHeatmap's pheatmap translation function.
}}

\section{Version 1.1.10}{\itemize{
\item Added rasterization for `dittoPlot()` and `dittoPlotVarsAcrossGroups()` under control by `do.raster` and `raster.dpi` inputs. (Rasterization effects jitter layers only.)
}}

\section{Version 1.1.9}{\itemize{
\item Added rasterization for `dittoScatterPlot()` and `dittoDimPlot()` under control by `do.raster` and `raster.dpi` inputs. Rasterization flattens the dots layer which can make working with a plot that has many, many cells, within a vector-based graphics editor, easier.
}}

\section{Version 1.1.8}{\itemize{
\item Added `dittoDotPlot()` function for suscinct features x groups visualization.
\item Expanded test coverage for multi-plotters.
\item Sped up `dittoPlotVarsAcrossGroups()`
\item Deprecated input: `OUT.List` was renamed to `list.out` in all instances (Affects: `multi_dittoPlot()`, `multi_dittoDimPlot()`, and `multi_dittoDimPlotVaryCells()`).
}}

\section{Version 1.1.7}{\itemize{
\item Added `order` input to dittoDimPlot and dittoScatterPlot for giving control over plotting order based on the colors data.
}}

\section{Version 1.1.6}{\itemize{
\item Updated the vignette to include Hex plotters and improve clarity.
\item Added `do.label` & `do.ellipse` features to Hex plotters.
}}

\section{Version 1.1.5}{\itemize{
\item Added features available in `dittoDimPlot()/Hex()` to `dittoScatterPlot()/Hex()`, including: `do.label`, `do.letter`, `do.ellipse`, `add.trajectory.lineages`, and `add.trajectory.curves`.
}}

\section{Version 1.1.4}{\itemize{
\item Bug fix: Allow `do.label` to work in conjunction with faceting of dittoDimPlots. Added a `labels.split.by` input in order to allow proper handling of group center calculations even when users manually add faceting via `extra.vars` instead of `split.by`.
}}

\section{Version 1.1.3}{\itemize{
\item Bug fix: Allow expression grabbing across visualizations (ultimately via `gene()`) from non-default assays of Seurat objects.
\item Added ignorance of `do.label`, `do.letter`, and `do.ellipse` inputs of dittoDimPlot for continuous target data, which allows their use in `multi_dittoDimPlot()` when `vars` includes both discrete & continuous data.
}}

\section{Version 1.1.2}{\itemize{
\item Added density plotting via hex binning with 2 new functions `dittoDimhex()` and `dittoScatterHex()`.
\item Added density display through contour line display for both the new "`Hex()`" plotters and their "`Plot()`" counterparts. This functionality is controlled via `do.contour`, `contour.color`, and `contour.linetype` inputs.
}}

\section{Version 1.1.1}{\itemize{
\item Added ability to adjust metadata during retrieval with `meta()`
\item Added ability to display metadata in `dittoHeatmaps()` (in addition to, or instead of, showing genes)
}}

\section{Version 1.0.0}{\itemize{
\item Submitted to Bioconductor.
}}
