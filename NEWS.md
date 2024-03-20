# dittoSeq 1.15.7

* Added new 'color.method = "prop.<value>"' functionality to dittoScatterHex() and dittoDimHex() which allows coloring bins by the proportion of a given 'color.var'-data value.

# dittoSeq 1.15.6

* Added ability to have 'dittoDotPlot()' group shown features ("vars") into categories, by providing 'vars' as a named list.  Also adds inputs 'categories.split.adjust' and 'categories.theme.adjust' which can be used to turn off the automatic addition of display-style adjustments that make category labels appear more like category labels.
* Added native control to 'dittoDotPlot()' over which axis is used for 'vars' versus 'group.by'-groupings, controlled by the new input: 'vars.dir'.
* BugFix: Corrects a bug in how the version of the Seurat-package was checked during expression matrix retrieval.

# dittoSeq 1.15.5

* BugFix: Corrects a bug which blocked provision of 'annotation_row' and 'annotation_colors' inputs to dittoHeatmap without also generating column annotations via either 'annot.by' or 'annotation_col'.

# dittoSeq 1.15.4

* Various behind the scenes updates to make warnings from ggplot and Seurat go away

# dittoSeq 1.15.3

* Added 'labels.repel.adjust' input to 'dittoDimPlot()', 'dittoScatterPlot()', 'dittoDimHex()', and 'dittoScatterHex()'. The input allows extra control of how labels will be placed when 'do.label = TRUE' (and 'labels.repel' is not set to FALSE) by providing a mechanism to pass desired parameters through to the ggrepel function used for plotting the labels.

# dittoSeq 1.15.2

* Added multi-modality feature expression retrieval to provide for across modality plotting. The 'assay' input can now be provided a vector of assays to target. The 'swap.rownames' input was also updated for multi-modality access purposes. The full system is documented in its own page, '?GeneTargeting', that function documentation for 'assay', 'slot', and 'swap.rownames' inputs have been updated to point to.
* Completed the deprecation of 'dittoHeatmap()'s 'highlight.genes' input

# dittoSeq 1.15.1

* Added 'vlnplot.quantiles' and 'boxplot.outlier.size' inputs to 'dittoPlot()', 'dittoPlotVarsAcrossGroups()', and 'dittoFreqPlot()' functions.

# dittoSeq 1.14

* Feature Extensions:
  1. 'dittoDotPlot()' & 'dittoPlotVarsAcrossGroups()': Improved 'group.by' ordering control via retention of factor levels and addition of a new 'groupings.drop.unused' input to control retention of empty levels.
  2. 'dittoHeatmap()': Targeting Seurat clusters with the "ident" shortcut now works for the 'annot.by' input of 'dittoHeatmap()'.
* Bug Fixes:
  1. 'dittoHeatmap()': Fixed a feature exclusion check in 'dittoHeatmap()' meant to remove features without any non-zero values. Previously, it removed all features with a mean of zero, which blocked plotting from pre-scaled data.
  2. 'dittoDimPlot()' & 'getReductions()': Eliminated cases where 'getReductions()' did not return NULL for 'object's containing zero dimensionality reductions. This fix also improves associated error messaging of 'dittoDimPlot()' where such cases were missed.

# dittoSeq 1.12

* No code updates. (Bioconductor version number updates only)

# dittoSeq 1.10

* Added ability to plot multiple 'var' in a single 'dittoPlot()', 'dittoDimPlot()', 'dittoScatterPlot()', 'dittoDimHex()', and 'dittoScatterHex()' call by giving a vector of genes or continuous metadata to the 'var' or 'color.var' input. Customization of how the "multivar" data is displayed can be controlled with:
1- 'multivar.aes' (context: 'dittoPlot()' only) - which plot aesthetic is utilized for displaying var-values.
2- 'multivar.split.dir' - faceting direction to use for var-data when combining with an additional 'split.by' variable.
* Improved the compatibility with 'split.by'/faceting customizations, specifically with 'split.adjust = list(scales = "free")', by making implementations of 'min'/'max' inputs less intrusive. Note: This change very minorly alters the default output of some plotters.
* Improved error messaging for cases where 'object' does not have cell/column names.

# dittoSeq 1.8

* Minor Feature Add: 'randomize' option for 'order' input of 'dittoDimPlot()' and 'dittoScatterPlot()'

# dittoSeq 1.6

* Vignette Update: Added a 'Quick-Reference: Seurat<=>dittoSeq' section.
* Build & Test Infrastructure Update: Removed Seurat dependency from all build and test materials by removing Seurat code from the vignette and making all unit-testing of Seurat interactions conditional on both presence of Seurat and successful SCE to Seurat conversion.
* Bug Fixes:
1- Fixed dittoFreqPlot calculation machinery to properly target all cell types but only necessary groupings for every sample. Removed the 'retain.factor.levels' input because proper calculations treat 'var'-data as a factor, and groupings data as non-factor.
2- Allowed dittoHeatmap() to properly 'drop_levels' of annotations by ensuring 'annotation_colors' is not populated with colors for empty levels which would be dropped.
3- Made 'do.label' machinery of scatter plots robust to NAs.

# dittoSeq 1.4

* Added 1 new Visualization function: 'dittoFreqPlot()'.
* Added interaction with 'rowData' of SE and SCEs via a 'swap.rownames' input, e.g. to simplify provision of 'var's via symbols vs IDs.
* Improved & expanded 'split.by' capabilities by:
1- adding them to 'dittoBarPlot()', 'dittoDotPlot()', and 'dittoPlotVarsAcrossGroups()';
2- adding 'split.adjust' input to all functions for passing adjustments to underlying 'facet_grid()' and 'facet_wrap()' calls;
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
