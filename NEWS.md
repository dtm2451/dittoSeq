# dittoSeq 1.8

* Minor Feature Add: 'randomize' option for 'order' input of 'dittoDimPlot()' and 'dittoScatterPlot()'

# dittoSeq 1.6

* Vignette Update: Added a 'Quick-Reference: Seurat<=>dittoSeq' section.
* Build & Test Infrastructure Update: Removed Seurat dependency from all build and test materials by removing Seurat code from the vignette and making all unit-testing of Seurat interactions conditional on both presence of Seurat and successful SCE to Seurat cnversion.
* Bug Fixes:
1- Fixed dittoFreqPlot calculation machinery to properly target all cell types but only necessary groupings for every sample. Removed the 'retain.factor.levels' input because proper calculations treat 'var'-data as a factor, and groupings data as non-factor.
2- Allowed dittoHeatmap() to properly 'drop_levels' of annotations by ensuring 'annotation_colors' is not populated with colors for empty levels which would be dropped.
3- Made 'do.label' machinery of scatter plots robust to NAs.

# dittoSeq 1.4

* Added 1 new Visualization function: 'dittoFreqPlot()'.
* Added interaction with 'rowData' of SE and SCEs via a 'swap.rownames' input, e.g. to simplify provision of 'var's via symbols vs IDs.
* Improved & expanded 'split.by' capabilities by:
1- adding them to 'dittoBarPlot()', 'dittoDotPlot()', and 'dittoPlotVarsAcrossGroups()';
2- adding 'split.adjust' input to all functions for passing adjudstments to underlying 'facet_grid()' and 'facet_wrap()' calls;
3- adding 'split.show.all.others' input to 'dittoDimPlot()' and 'dittoScatterPlot()' to allow the full spectrum of points, rather than just points excluded with 'cells.use', to be shown as light gray in the background of all facets;
4- Bug fix: splitting now works with labeling of Dim/Scatter plots, with label position calculated per facet, and without affecting facet order.
* Improved 'dittoPlot()'-plotting engine (also effects 'dittoPlotVarsAcrossGroups()', and 'dittoFreqPlot()') by:
for y-axis plotting,
1- extended geom dodging to also work on jitters when 'color.by' is used to add subgroupings &
2- added a 'boxplot.lineweight' control option;
for x-axis / ridge plotting,
1- added an alternative histogram-shaping option (Try 'ridgeplot.shape = "hist"') &
2- improved use of white space via a new 'ridgeplot.ymax.expansion' input.
* Standardized output logic so that 'do.hover = TRUE' will lead to plotly conversion even when 'data.out = TRUE'. 
* 'dittoHeatmap()': 'order.by' can also now accept multiple gene/metadata names to order by & bug fix: when given an integer vector, that vector will be used directly to set the order of heatmap columns.
* 'dittoBarPlot()': grouping & 'var' order control improved via addition of a 'retain.factor.levels' input.

# dittoSeq 1.2

* Added 3 New Visualization Functions, 'dittoDotPlot()', 'dittoDimHex()' & 'dittoScatterHex()'.
* Expanded SummarizedExperiment compatibility across the entire toolset.
* Added ComplexHeatmap integration to 'dittoHeatmap()', controlled by a new input, 'complex'.
* Added Rasterization for improved image editor compatibility of complex plots. (See the dedicated section in the vignette for details.)
* Added 'labels.split.by' input & 'do.contour', 'contour.color', and 'contour.linetype' inputs to scatter/dim-plots.
* Added 'order' input to scatter/dim-plots for control of plotting order.
* Added 'metas' input for displaying such data with 'dittoHeatmap()'.
* Added 'adjustment' input to 'meta()', which works exactly as in 'gene()' (but this is not yet implemented within data grab of visualiation functions).
* Added 'adj.fxn' input to 'meta()' aand 'gene()' for added control of how data might be adjusted (but this is not yet implemented within data grab of visualiation functions).
* Replaced (deprecated) 'highlight.genes' input with 'highlight.features' in 'dittoHeatmap()'.
* Replaced (deprecated) 'OUT.List' input with 'list.out' for all 'multi_*' plotters.


# dittoSeq 1.0

* Submitted to Bioconductor.
