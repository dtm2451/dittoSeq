# dittoSeq 1.17.1

* Minor changes to ensure compatibility with ggplot2 version 4.0.0 changes. Notably, the 'MASS' package is now added as a suggested companion package, and error messages suggesting its installation were added for features which rely on functionality from this package. This change is concurrent with ggplot2 itself downgrading the tool from an imported requirement to a suggested addition.

# dittoSeq 1.16

* Feature Extensions:
    1. Multi-modality functionality: To support visualization of markers from multiple modalities in the same plot, e.g. gene expression by RNA and protein capture by ADT, the mechanics of 'assay', 'slot', and 'swap.rownames' inputs have been expanded, although defaults are unchanged. See the '?GeneTargeting' documentation page for details. For the standard Seurat CITEseq case, set 'assay = c("RNA", "ADT")'.
    2. 'dittoDotPlot()' vars-categories: Added support for categorization of markers, as well as for x and y axes swapping.
        * Provide 'vars' as a named list to group markers (list element values) into categories (list element names).
        * New input 'vars.dir' controls which axis is used for markers ("x" by default, or "y").
        * New boolean inputs 'categories.theme.adjust' or 'categories.split.adjust' can be used to turn off associated automated additions to the 'theme' input or 'split.adjust' input as well as faceting mechanics, respectively.
    3. 'dittoDotPlot()' 3-color scaling: Added support for injecting a midpoint color to the color scale via 2 new inputs.
        * New input 'mid.color' acts as the switch, and can be set to the specific strings "ryb", "rwb", or "rgb" (g for gray here) for a single-point quick update to use corresponding ColorBrewer inspired scales (effectively updating 'min.color' and 'max.color' as well). 'mid.color' can alternatively be given a color directly for more fine-grain control of colors.
        * New input 'mid' controls the data value at which 'mid.color' will be used in the scale.
    4. Additional data representation controls for all 'dittoPlot()'-style plotters, which includes 'dittoFreqPlot()':
        * New input 'boxplot.outlier.size' allows control of the outlier shape's size for "boxplot" representations.
        * New input 'vlnplot.quantiles' allows addition of lines at requested data quantiles for "vlnplot" representations.
    5. Added a new built in 'color.method' style for 'dittoScatterHex()' and 'dittoDimHex()' plotters:
        * When 'color.var' targets discrete data, giving 'color.method = "prop.\<value\>"', where \<value\> is an actual data level of 'color.var'-data, will set coloring to represent the proportion of \<value\> among the 'color.var'-data of each bin.
    6. New input 'labels.repel.adjust' allows finer control of the 'do.label' plot addition, via input pass-through to the geom functions underlying labeling. This affects 'dittoDimPlot()', 'dittoScatterPlot()', 'dittoDimHex()', and 'dittoScatterHex() functions.
* Bug Fixes:
    1. 'dittoHeatmap()': Fixed a bug which blocked provision of 'annotation_row' and 'annotation_colors' inputs to 'dittoHeatmap()' without also generating column annotations via either 'annot.by' or direct 'annotation_col' provision.
* Deprecation:
    1. Completed deprecation of 'dittoHeatmap()'s 'highlight.genes' input via removal from the function.
* Dependency Upkeep (generally invisible to users):
    1. ggplot-v3: Replaced all calls to the deprecated 'aes_string()' function with calls to the standard 'aes()' function. In cases where mappings are successively built internally to accommodate customization or flexibility, 'modifyList()' usage replaces the previous simple 'list' and 'do.call()' management.
    2. Seurat-v5: When the user's Seurat package version is 5.0 or higher, conditional code switches expression data retrieval from a call to the reportedly superseded 'GetAssayData()' function to the newly supported 'SeuratObj[[\<assay\>]][\<slot\>]' syntax.

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
