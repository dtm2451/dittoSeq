\name{dittoSeqNews}
\title{dittoSeq News}
\encoding{UTF-8}

\section{Version 1.3.13}{\itemize{
\item Adds a `split.adjust` input, to all functions with `split.by`, for passing adjudstments like 'scales = "free"' through in the code where faceting is added.
\item Documentation bug fix
}}

\section{Version 1.3.12}{\itemize{
\item Adds the entirely new Viz fxn: `dittoFreqPlot()`.
\item Improves the `dittoPlot()`-plotting engine (which affects `dittoPlot()` & `dittoPlotVarsAcrossGroups()`, and `dittoFreqPlot()`) by: incorporating a position-dodge fix for jitter-plots when `color.by` is used for subgrouping; and adding a `boxplot.lineweight` control option.
}}

\section{Version 1.3.11}{\itemize{
\item Adds `split.by` faceting to functions which did not have it: `dittoBarPlot()`, `dittoDotPlot()`, and `dittoPlotVarsAcrossGroups()`.
\item Improves `split.by` for `dittoDimPlot()` and `dittoScatterPlot()` by adding a `split.show.all.others` input. The input controls whether the full spectrum of points, versus just points excluded with `cells.use`, will be shown as light gray in the background of plots facets.
}}

\section{Version 1.3.10}{\itemize{
\item Standardizes `data.out` & `do.hover` interplay, for all visualizations with `do.hover`, to allow both plotly conversion and data output.
}}

\section{Version 1.3.9}{\itemize{
\item Adds an alternative histogram-shaping to ridge-plots via 3 new inputs: `ridgeplot.shape = "hist"` to turn it on, and `ridgeplot.bins` and `ridgeplot.binwidth` to adjust the look.
}}

\section{Version 1.3.8}{\itemize{
\item Build-side Bug Fix, minor: Added `rmarkdown` to the list of `Suggests:` packages.
}}

\section{Version 1.3.7}{\itemize{
\item Bug Fix, minor: Remove link to Seurat's `pbmc_small` object which moved to their SeuratObject package.
}}

\section{Version 1.3.6}{\itemize{
\item Bug Fix: numeric vectors given to `dittoHeatmap()` `order.by` input are now used directly to set order of heatmap columns, as originally intended.
\item Feature Update: `dittoHeatmap()`'s `order.by` can also now, alternatively, accept multiple gene/metadata names, rather than just one.
}}

\section{Version 1.3.5}{\itemize{
\item Unit test correction
}}

\section{Version 1.3.4}{\itemize{
\item Adds citation for the Advance Access version of the dittoSeq Applications Note published with Bioinformatics.
\item Updates the `do.hover` logic for dittoPlot & dittoPlotVarsAcrossGroups so that plotly conversion is skipped for ridgeplots, yet otherwise is now performed even when `data.out = TRUE` or when no jitter is requested.
}}

\section{Version 1.3.3}{\itemize{
\item Adds a `swap.rownames` input, to allow utilization of SE/SCE rowData values, in all visualizations that can target feature data, as well as in the `gene()`, `isGene()`, and `getGenes()` helper functions.
\item Update the `importDittoBulk()` function to remove the edgeR package dependency.
}}

\section{Version 1.3.2}{\itemize{
\item Bug Fix: keep `split.by` factor level ordering when used with `do.label`
}}

\section{Version 1.3.1}{\itemize{
\item Bug Fix: Better y-axis expansion for `dittoPlots()` with many ridgeplots
}}

\section{Version 1.2}{\itemize{
\item Added 3 New Visualization Functions, `dittoDotPlot()`, `dittoDimHex()` & `dittoScatterHex()`.
\item Expanded SummarizedExperiment compatibility across the entire toolset.
\item Added ComplexHeatmap integration to `dittoHeatmap()`, controlled by a new input, `complex`.
\item Added Rasterization for improved image editor compatibility of complex plots. (See the dedicated section in the vignette for details.)
\item Added `labels.split.by` input & `do.contour`, `contour.color`, and `contour.linetype` inputs to scatter/dim-plots.
\item Added `order` input to scatter/dim-plots for control of plotting order.
\item Added `metas` input for displaying such data with `dittoHeatmap()`.
\item Added `adjustment` input to `meta()`, which works exactly as in `gene()` (but this is not yet implemented within data grab of visualiation functions).
\item Added `adj.fxn` input to `meta()` aand `gene()` for added control of how data might be adjusted (but this is not yet implemented within data grab of visualiation functions).
\item Replaced (deprecated) `highlight.genes` input with `highlight.features` in `dittoHeatmap()`.
\item Replaced (deprecated) `OUT.List` input with `list.out` for all `multi_*` plotters.
}}

\section{Version 1.0.0}{\itemize{
\item Submitted to Bioconductor.
}}
