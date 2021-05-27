\name{dittoSeqNews}
\title{dittoSeq News}
\encoding{UTF-8}

\section{Version 1.5.1}{\itemize{
\item Build & Test Infrastructure Fix/Update: Removes the remaining Seurat code from the vignette and makes all unit-testing of Seurat interactions conditional on SCE to Seurat cnversion.
\item Vignette Update: Adds a 'Quick-Reference: Seurat<=>dittoSeq' section
\item Doc Update: 'metaLevels()'
}}

\section{Version 1.4}{\itemize{
\item Added 1 new Visualization function: 'dittoFreqPlot()'.
\item Added interaction with 'rowData' of SE and SCEs via a 'swap.rownames' input, e.g. to simplify provision of 'var's via symbols vs IDs.
\item Improved & expanded 'split.by' capabilities by:
1- adding them to 'dittoBarPlot()', 'dittoDotPlot()', and 'dittoPlotVarsAcrossGroups()';
2- adding 'split.adjust' input to all functions for passing adjudstments to underlying 'facet_grid()' and 'facet_wrap()' calls;
3- adding 'split.show.all.others' input to 'dittoDimPlot()' and 'dittoScatterPlot()' to allow the full spectrum of points, rather than just points excluded with 'cells.use', to be shown as light gray in the background of all facets;
4- Bug fix: splitting now works with labeling of Dim/Scatter plots, with label position calculated per facet, and without affecting facet order.
\item Improved 'dittoPlot()'-plotting engine (also effects 'dittoPlotVarsAcrossGroups()', and 'dittoFreqPlot()') by:
for y-axis plotting,
1- extended geom dodging to also work on jitters when 'color.by' is used to add subgroupings &
2- added a 'boxplot.lineweight' control option;
for x-axis / ridge plotting,
1- added an alternative histogram-shaping option (Try 'ridgeplot.shape = "hist"') &
2- improved use of white space via a new 'ridgeplot.ymax.expansion' input.
\item Standardized output logic so that 'do.hover = TRUE' will lead to plotly conversion even when 'data.out = TRUE'. 
\item 'dittoHeatmap()': 'order.by' can also now accept multiple gene/metadata names to order by & bug fix: when given an integer vector, that vector will be used directly to set the order of heatmap columns.
\item 'dittoBarPlot()': grouping & 'var' order control improved via addition of a 'retain.factor.levels' input.
}}

\section{Version 1.2}{\itemize{
\item Added 3 New Visualization Functions, 'dittoDotPlot()', 'dittoDimHex()' & 'dittoScatterHex()'.
\item Expanded SummarizedExperiment compatibility across the entire toolset.
\item Added ComplexHeatmap integration to 'dittoHeatmap()', controlled by a new input, 'complex'.
\item Added Rasterization for improved image editor compatibility of complex plots. (See the dedicated section in the vignette for details.)
\item Added 'labels.split.by' input & 'do.contour', 'contour.color', and 'contour.linetype' inputs to scatter/dim-plots.
\item Added 'order' input to scatter/dim-plots for control of plotting order.
\item Added 'metas' input for displaying such data with 'dittoHeatmap()'.
\item Added 'adjustment' input to 'meta()', which works exactly as in 'gene()' (but this is not yet implemented within data grab of visualiation functions).
\item Added 'adj.fxn' input to 'meta()' aand 'gene()' for added control of how data might be adjusted (but this is not yet implemented within data grab of visualiation functions).
\item Replaced (deprecated) 'highlight.genes' input with 'highlight.features' in 'dittoHeatmap()'.
\item Replaced (deprecated) 'OUT.List' input with 'list.out' for all 'multi_*' plotters.
}}

\section{Version 1.0.0}{\itemize{
\item Submitted to Bioconductor.
}}
