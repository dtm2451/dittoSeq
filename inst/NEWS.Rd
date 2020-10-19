\name{dittoSeqNews}
\title{dittoSeq News}
\encoding{UTF-8}

\section{Version 1.2}{\itemize{
\item Added 3 New Visualization Functions, `dittoDotPlot()`, `dittoDimScatter()` & `dittoDimHex()`.
\item Expanded SummarizedExperiment compatibility across the entire toolset.
\item Added ComplexHeatmap integration to `dittoHeatmap()`, controlled by a new input, `complex`.
\item Added Rasterization for improved image editor compatibility of complex plots. (See the dedicated section in the vignette for details.)
\item Added `labels.split.by` input & `do.contour`, `contour.color`, and `contour.linetype` inputs to scatter/dim-plots.
\item Added `order` input to scatter/dim-plots for control of plotting order.
\item Added `metas` input for displaying such data with `dittoHeatmap()`.
\item Replaced (deprecated) `highlight.genes` input with `highlight.features` in `dittoHeatmap()`.
\item Replaced (deprecated) `OUT.List` input with `list.out` for all `multi_*` plotters.
}}

\section{Version 1.0.0}{\itemize{
\item Submitted to Bioconductor.
}}
