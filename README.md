# <img src="vignettes/dittoSeq_HexSticker.png" alt="dittoSeq" height="200"> dittoSeq 

**A set of functions built to enable analysis and visualization of single-cell and bulk RNA-sequencing data by novice, experienced, and color blind coders**

dittoSeq includes universal plotting and helper functions for working with (sc)RNAseq data processed in these packages:

- single-cell:
  - Seurat (versions 2+), *Seurat* data structure
  - scran / scater / other Bioconductor packages that utilize the *SingleCellExperiment* data structure
- bulk:
  - edgeR, *DGEList* data structure
  - DESeq2 / other Bioconductor packages that utilize the *SummarizedExperiment* data structure
  
All plotting functions spit out easy-to-read, color blind friendly, plots (ggplot2, plotly, or pheatmap/ComplexHeatmap) upon minimal coding input for your daily analysis needs, yet also allow sufficient manipulations to provide for out-of-the-box submission-quality figures!

dittoSeq also makes access of underlying data easy, for submitting to journals or for adding extra layers to the plot, with `data.out = TRUE` inputs!

![Overview](vignettes/dittoSeq.gif)

### News:

**Major functionality updates are coming in the next release!**

Current updates in 'devel' dittoSeq v1.15, coming soon in Bioconductor 3.19:

- Feature Extensions:

    1. **Multi-modality functionality**: To support visualization of markers from multiple modalities within a single plot, the way that `assay` and `slot` inputs can be provided for Seurat-objects, and that `assay` and `swap.rownames` can be provided for SingleCellExperiment and Summarized Experiment objects has been overhauled.  A new documentation page, [`?GeneTargeting`](man/GeneTargeting.Rd), describes the new methodologies.
        - Note: *Single assay mode currently remains the default for all plotters*. In the current implementation, you need to explicitly set `assay` whenever aiming to target multiple modalities, but I will be considering defaulting to e.g. `assay = c("RNA", "ADT")` for Seurat objects in a future dittoSeq-v2.0.
    2. **'dittoDotPlot()' vars-categories**: Added support for categorization of `vars`/markers shown in 'dittoDotPlot()'s and also for axes swapping.
        - Added `vars.dir` input to give internal control over whether markers are shown on the x-axis (default) or the y-axis (`vars.dir = "y"`).
        - `vars` input can be given as a named list to group markers (list element values) into categories (list element names).
        - Added automatic addition of display-style adjustments that make category labels appear more like category labels.
            - New inputs `categories.split.adjust` and `categories.theme.adjust` were added to let users turn off these display adjustments. Elements which are added to `split.adjust` (and then ultimately given to `ggplot2::facet_grid()`) can be turned off by setting `categories.split.adjust = FALSE` and elements which are added to `theme` (and applied via `ggplot2::theme()`) can be turned off by setting `categories.theme.adjust = FALSE`.
    3. (Underway in #147) **'dittoDotPlot()' 3-color scale**: Added support for adding a midpoint color to the color scale used for 'dittoDotPlot()'s.
        - New input `mid.color` controls the switch:
            - Left as the default, `NULL`, a 2-color scale is used (from 'min.color' to 'max.color')
            - Given `mid.color = "<color>"`, a 3-color scale is used (from 'min.color' to \<color\> to 'max.color')
            - Giving `mid.color = "rgb"` or `"rwb"` allows single-point quick update to all of 'min.color', 'mid.color', and 'max.color' for use of one of two standard 3-color scales ("rgb": from "blue" to "gray90" to "red"; "rwb": from "blue" to "white" to "red").
        - New input `mid` controls the data value at which 'mid.color' will be used in the scale, and receives intuitive defaulting so users generally don't need to provide it.
        - *This mechanism is being rolled out an tested with `dittoDotPlot()` first, but* ***users can expect extension of this functionality to other visualizations in an upcoming release!***
    4. 'dittoDimPlot()', 'dittoScatterPlot()', 'dittoDimHex()', and 'dittoScatterHex():
        - Added 'labels.repel.adjust' input which provide additional control of labeling via input pass-through to the `ggrepel::geom_label_repel()` (`labels.highlight = TRUE`, the default) and `ggrepel::geom_text_repel()` (`labels.highlight = FALSE`) functions which underly `do.label = TRUE` labeling (when `labels.repel` is left as the default `TRUE`)
    5. 'dittoPlot()', 'dittoFreqPlot()', and 'dittoPlotVarsAcrossGroups()':
        - Added 'boxplot.outlier.size' input to allow control of the outlier shape's size
        - Added 'vlnplot.quantiles' input to allow drawing of lines, within violin plot data representations, at requested data quantiles.
    6. 'dittoScatterHex()' and 'dittoDimHex()':
        - Added a new 'color.method' style for discrete 'color.var'-data. Users can give `color.method = "prop.<value>"`, where \<value\> is an actual data level of 'color.var'-data, to have color represent the proportion of `'color.var'-data == <value>` for all bins.

- Bug Fixes:
    
    1. 'dittoHeatmap()': Fixed a bug which blocked provision of 'annotation_row' and 'annotation_colors' inputs to dittoHeatmap without also generating column annotations via either 'annot.by' or direct 'annotation_col' provision. 

- Upkeep with ggplot-v3 & Seurat-v5, *details here are generally invisbile to the user*:
  
    1. ggplot-v3: Switched from making use of `do.call()` with the deprecated `aes_string()`, and simple `list` management for successively built setups, to mostly direct `aes(.data[[<col>]])` calls, and use of `modifyList` for additions in successively built setups.  This methodology should be backwards compatible to earlier ggplot versions, but that has not been officially tested.
    2. Seurat-v5: Added conditional code that switched to the newly supported `SeuratObj[[<assay>]][<slot>]` syntax for expression data retrieval when the user's Seurat package version is 5.0 or higher. 

#### Previous updates:

<details>

  <summary>Updates in dittoSeq v1.14 (Bioconductor 3.18)</summary>
  
  - Feature Extensions:
    
      1. 'dittoDotPlot()' & 'dittoPlotVarsAcrossGroups()': Improved 'group.by' ordering control via retention of factor levels and addition of a new 'groupings.drop.unused' input to control retention of empty levels.
      2. 'dittoHeatmap()': Targeting Seurat clusters with the "ident" shortcut now works for the 'annot.by' input of 'dittoHeatmap()'.
  
  - Bug Fixes:
    
      1. 'dittoHeatmap()': Fixed a feature exclusion check in 'dittoHeatmap()' meant to remove features without any non-zero values. Previously, it removed all features with a mean of zero, which blocked plotting from pre-scaled data.
      2. 'dittoHeatmap()': (New in dittoSeq-v1.15-devel, but also pushed to the released v1.14.2) Fixed a bug which blocked provision of 'annotation_row' and 'annotation_colors' inputs to dittoHeatmap without also generating column annotations via either 'annot.by' or direct 'annotation_col' provision.
      3. 'dittoDimPlot()' & 'getReductions()': Eliminated cases where 'getReductions()' did not return NULL for 'object's containing zero dimensionality reductions. This fix also improves associated error messaging of 'dittoDimPlot()' where such cases were missed.
  
  - Upkeep with ggplot-v3 & Seurat-v5, *details here are generally invisbile to the user* (New in dittoSeq-v1.15-devel, but also pushed to the released v1.14.1):
  
      1. ggplot-v3: Switched from making use of `do.call()` with the deprecated `aes_string()`, and simple `list` management for successively built setups, to mostly direct `aes(.data[[<col>]])` calls, and use of `modifyList` for additions in suucessively built setups.  This methodology should be backwards compatible to earlier ggplot versions, but that has not been officially tested.
      2. Seurat-v5: Added conditional code that switched to the newly supported `SeuratObj[[<assay>]][<slot>]` syntax for expression data retrieval when the user's Seurat package version is 5.0 or higher.
  
</details>


<details>

  <summary>No code updates in dittoSeq v1.10 & v1.12 (Bioconductor 3.16 & 3.17)</summary>
  
  - Bioconductor-maintained version number updates only
  
</details>


<details>

  <summary>Updates in dittoSeq v1.8 (Bioconductor 3.15)</summary>
  
  - Added 'randomize' option for 'order' input of 'dittoDimPlot()' and 'dittoScatterPlot()'
  
</details>


<details>

  <summary>Updates in dittoSeq v1.6 (Bioconductor 3.14)</summary>
  
  - Vignette Update: Added a 'Quick-Reference: Seurat<=>dittoSeq' section.
  - Build & Test Infrastructure Update: Removed Seurat dependency from all build and test materials by removing Seurat code from the vignette and making all unit-testing of Seurat interactions conditional on both presence of Seurat and successful SCE to Seurat cnversion.
  - Bug Fixes:

    1. Fixed dittoFreqPlot calculation machinery to properly target all cell types but only necessary groupings for every sample. Removed the 'retain.factor.levels' input because proper calculations treat 'var'-data as a factor, and groupings data as non-factor.
    2. Allowed dittoHeatmap() to properly 'drop_levels' of annotations by ensuring 'annotation_colors' is not populated with colors for empty levels which would be dropped.
    3. Made 'do.label' machinery of scatter plots robust to NAs.
  
</details>


<details>

  <summary>Updates in dittoSeq v1.4 (Bioconductor 3.13)</summary>

  - Added 1 New Visualization Function: `dittoFreqPlot()`:
    - Combines the population frequency summarization of `dittoBarPlot()` with the plotting style of `dittoPlot()` to enable per-population, per-sample, per-group frequency comparisons which focus on individual cell types / clusters!
  - Improved & expanded faceting capabilities with `split.by` inputs:
      - Added `split.by` to functions which did not have it: `dittoBarPlot()`, `dittoDotPlot()`, and `dittoPlotVarsAcrossGroups()` 
      - Added `split.adjust` input to allow tweaks to the underlying `facet_grid()` and `facet_wrap()` calls.
      - Better compatibility with other features
          - works with labeling of Dim/Scatter plots
          - new `split.show.all.others` input now controls whether the full spectrum of points, versus just points excluded with `cells.use`, will be shown as light gray in the background of Dim/Scatter facets.
  - Improved `dittoPlot()`-plotting engine:
      - y-axis plotting:
          - geom dodging when `color.by` is used to add subgroupings now works for jitters too.
          - added a `boxplot.lineweight` control option.
      - x-axis / ridge-plotting:
          - Added an alternative histogram-shaping option (Try adding `ridgeplot.shape = "hist"`!)
          - Better use of white space (via adjustments to default plot grid expansion & exposure of a `ridgeplot.ymax.expansion` input to allow user override.)
  - Improved ordering capability for `dittoHeatmap()` & `dittoBarPlot()`:
      - `dittoHeatmap()`: You can now give many metadata to `order.by` and it will use them all, prioritizing earliest items
      - `dittoBarPlot()`: Factor-level ordering can now be retained in dittoBarPlot for `var` and `group.by` data, a typically expected behavior, by setting a new input `retain.factor.levels = TRUE`.
  - Added interaction with `rowData` of SE and SCEs:
      - `swap.rownames` input allows indication of genes/rows by non-default rownames. E.g. for an `object` with Ensembl_IDs as the default and a rowData column named 'symbol' that contains gene symbols, those symbols can be used via `dittoFunction(..., var = "<gene_symbol>", swap.rownames = "symbol"`).
  - Quality of Life improvements:
      - Standardized `data.out` & `do.hover` interplay to allow both plotly conversion and data output.
      - Documentation Updates
  
</details>


<details>

  <summary>Updates in dittoSeq v1.2 (Bioconductor 3.12)</summary>
  
  - Added 3 New Visualization Functions, `dittoDotPlot()`, `dittoDimHex()` & `dittoScatterHex()`.
  - Expanded SummarizedExperiment compatibility across the entire toolset.
  - Added ComplexHeatmap integration to `dittoHeatmap()`, controlled by a new input, `complex`.
  - Added Rasterization for improved image editor compatibility of complex plots. (See the dedicated section in the vignette for details.)
  - Added `labels.split.by` input & `do.contour`, `contour.color`, and `contour.linetype` inputs to scatter/dim-plots.
  - Added `order` input to scatter/dim-plots for control of plotting order.
  - Added `metas` input for displaying such data with `dittoHeatmap()`.
  - Added `adjustment` input to `meta()`, which works exactly as in `gene()` (but this is not yet implemented within data grab of visualization functions).
  - Added `adj.fxn` input to `meta()` and `gene()` for added control of how data might be adjusted (but this is not yet implemented within data grab of visualization functions).
  - Replaced (deprecated) `highlight.genes` input with `highlight.features` in `dittoHeatmap()`.
  - Replaced (deprecated) `OUT.List` input with `list.out` for all `multi_*` plotters.
  
</details>


<details>

  <summary>Updates in dittoSeq v1.0 (Bioconductor 3.11)</summary>
  
  - Submitted to Bioconductor
  
</details>


### Color Blindness Compatibility:

The default colors of this package are meant to be color blind friendly.  To make it so, I used the suggested colors from this source: [Wong B, "Points of view: Color blindness." Nature Methods, 2011](https://www.nature.com/articles/nmeth.1618) and adapted them slightly by appending darker and lighter versions to create a 24 color vector. All plotting functions use these colors, stored in `dittoColors()`, by default. Also included is a Simulate() function that allows you to see what your function might look like to a colorblind individual. For more info on that, see the [Color blindness Friendliness section below](#color-blindness-friendliness)

### Demuxlet Tools

Included in this package currently are a set of functions to facilitate Mux-seq applications. For information about how to use these tools, see the [Demuxlet section down below](#demuxlet-tools). For more information on Demuxlet and Mux-sequencing, see the [Demuxlet GitHub Page](https://github.com/statgen/demuxlet). (Impetus: Many Mux-seq experiments will involve generating the side-by-side bulk and single-cell RNAseq data like the rest of the package is built for.)

## Installation:

```
### For R-4.0 users:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dittoSeq")

### For users with older versions of R:
# BiocManager will not let you install the pre-compiled version, but you can
# install directly from this GitHub via:
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("dtm2451/dittoSeq")
```

# Quick Reference: Seurat <=> dittoSeq

Because often users will be familiar with Seurat already, so this may be 90% of what you may need!

<details>

  <summary>Click to expand</summary>

  As of May 25th, 2021, Seurat-v4.0.2 & dittoSeq v1.4.1
  
  **Functions**
  
  Seurat Viz Function(s) | dittoSeq Equivalent(s)
  --- | ---
  DimPlot/ (I)FeaturePlot / UMAPPlot / etc. | dittoDimPlot / multi_dittoDimPlot
  VlnPlot / RidgePlot | dittoPlot / multi_dittoPlot
  DotPlot | dittoDotPlot
  FeatureScatter / GenePlot | dittoScatterPlot
  DoHeatmap | dittoHeatmap*
  [No Seurat Equivalent] | dittoBarPlot / dittoFreqPlot
  [No Seurat Equivalent] | dittoDimHex / dittoScatterHex
  [No Seurat Equivalent] | dittoPlotVarsAcrossGroups
  SpatialDimPlot, SpatialFeaturePlot, etc. | dittoSpatial (coming soon!)
  
  *Not all dittoSeq features exist in Seurat counterparts, and occasionally the
  same is true in the reverse.
  
  **Inputs**
  
  See reference below for the equivalent names of major inputs
  
  Seurat has had inconsistency in input names from version to version. dittoSeq
  drew some of its parameter names from previous Seurat-equivalents to ease
  cross-conversion, but continuing to blindly copy their parameter standards will
  break people's already existing code. Instead, dittoSeq input names are
  guaranteed to remain consistent across versions, unless a change is required for
  useful feature additions.
  
  Seurat Viz Input(s) | dittoSeq Equivalents
  --- | ---
  `object` | SAME
  `features` | `var` / `vars` (generally the 2nd input, so name not needed!) OR `genes` & `metas` for dittoHeatmap()
  `cells` (cell subsetting is not always available) | `cells.use` (consistently available)
  `reduction` & `dims` | `reduction.use` & `dim.1`, `dim.2`
  `pt.size` | `size` (or `jitter.size`)
  `group.by` | SAME
  `split.by` | SAME
  `shape.by` | SAME and also available in dittoPlot()
  `fill.by` | `color.by` (can be used to subset `group.by` further!)
  `assay` / `slot` | SAME
  `order` = logical | `order` but = "unordered" (default), "increasing", or "decreasing"
  `cols` | `color.panel` for discrete OR `min.color`, `max.color` for continuous
  `label` & `label.size` & `repel` | `do.label` & `labels.size` & `labels.repel`
  `interactive` | `do.hover` = via plotly conversion
  [Not in Seurat] | `data.out`, `do.raster`, `do.letter`, `do.ellipse`, `add.trajectory.lineages` and others!
  
</details>

# Quick Start Guide:

Load in your data, then go!:

```
library(dittoSeq)

# dittoSeq works natively with Seurat, SingleCellExperiment (SCE),
#   & SummarizedExperiment (SE) objects

# Seurat
seurat <- Seurat::pbmc_small
dittoPlot(seurat, "CD14", group.by = "ident")

# SingleCellEXperiment
sce <- Seurat::as.SingleCellExperiment(seurat)
dittoDimPlot(sce, "CD14")

# SummarizedExperiment
# (Please excuse the janky setup code for this quick example.)
library(SummarizedExperiment)
se <- as(as.SingleCellExperiment(Seurat::pbmc_small), "SummarizedExperiment")
rownames(se) <- rownames(sce)
dittoBarPlot(sce, "ident", group.by = "RNA_snn_res.0.8")

# For working with non-SE bulk RNAseq data, first import your data into a
#   SingleCellExperiment structure, (which is essentially a SummarizedExperiment
#   structure just with an added space for holding dimensionality reductions).
# myRNA <- importDittoBulk(dds) # DESeq2
# myRNA <- importDittoBulk(dgelist) # edgeR
# Then add dimensionality reductions
# myRNA <- addDimReduction(myRNA, embeddings, "pca")
#   above, embeddings = the dim-reduction matrix
myRNA <- example("importDittoBulk")

# You're ready!
dittoDimPlot("gene1", myRNA, size = 3)
```

Quickly determine the metadata and gene options for plotting with universal helper functions:

```
getMetas(seurat)
isMeta("nCount_RNA", seurat)

getGenes(myRNA)
isGene("CD3E", myRNA)

getReductions(sce)

# View them with these:
gene("CD3E", seurat, assay = "RNA", slot = "counts")
meta("groups", seurat)
metaLevels("groups", seurat)
```

### There are many dittoSeq Plot Types

**Intuitive default adjustments generally allow creation of immediately useable plots.**

```
# dittoDimPlot
dittoDimPlot(seurat, "ident", size = 3)
dittoDimPlot(seurat, "CD3E", size = 3)
```

![](vignettes/QuickStart1_Dim.png)

```
# dittoBarPlot
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8")
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    scale = "count")
```

![](vignettes/QuickStart2_Bar.png)

```
# dittoPlot
dittoPlot(seurat, "CD3E", group.by = "ident")
dittoPlot(seurat, "CD3E", group.by = "ident",
    plots = c("boxplot", "jitter"))
dittoPlot(seurat, "CD3E", group.by = "ident",
    plots = c("ridgeplot", "jitter"))
```

![](vignettes/QuickStart3_Plot.png)

```
# dittoHeatmap
dittoHeatmap(seurat, genes = getGenes(seurat)[1:20])
dittoHeatmap(seurat, genes = getGenes(seurat)[1:20],
    annot.by = c("groups", "nFeature_RNA"),
    scaled.to.max = TRUE,
    treeheight_row = 10)
# Turning off cell clustering can be necessary for large scRNAseq data
# Thus, clustering is turned off by default for single-cell data, but not for
# bulk RNAseq data.
# To control ordering/clustering separately, use 'order.by' or 'cluster_cols'
## (Not shown) ##
dittoHeatmap(seurat, genes = getGenes(seurat)[1:20],
    order.by = "groups")
dittoHeatmap(seurat, genes = getGenes(seurat)[1:20],
    cluster_cols = FALSE)
```

![](vignettes/QuickStart4_Heatmap.png)

```
# dittoScatterPlot
dittoScatterPlot(
    object = seurat,
    x.var = "CD3E", y.var = "nCount_RNA",
    color.var = "ident", shape.by = "RNA_snn_res.0.8",
    size = 3)
dittoScatterPlot(
    object = seurat,
    x.var = "nCount_RNA", y.var = "nFeature_RNA",
    color.var = "CD3E",
    size = 1.5)
```

![](vignettes/QuickStart5_Scatter.png)

```
# Also multi-plotters:
    # multi_dittoDimPlot (multiple, in an array)
    # multi_dittoDimPlotVaryCells (multiple, in an array, but showing only
    #     certain cells in each plot)
    # multi_dittoPlot (multiple, in an array)
    # dittoPlot_VarsAcrossGroups (multiple genes or metadata as the jitter
    #     points (and other representations), summarized across groups by
    #     z-score, or mean, or median, or any function that outputs a
    #     single numeric value from a numeric vector input.)
```

**Many adjustments can be made with simple additional inputs:**

dittoSeq allows many adjustments to how data is represented via inputs directly within dittoSeq functions.
Adjustments that are common across functions are briefly described below.
Some others are within the examples above.

For more details, review the full vignette (`vignette("dittoSeq")` after installation via Bioconductor)
and/or the documentation of individual functions (example: `?dittoDimPlot`).

Common Adjustments:

- All Titles are adjustable.
- Easily subset the cells shown with `cells.use`
- Colors can be adjusted easily.
- Underlying data can be output.
- plotly hovering can be added.
- Many more! (Legends removal, label rotation, labels' and groupings' names, ...)

```
# Adjust titles
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    main = "Starters",
    sub = "By Type",
    xlab = NULL,
    ylab = "Generation 1",
    x.labels = c("Ash", "Misty"),
    legend.title = "Types",
    var.labels.rename = c("Fire", "Water", "Grass"),
    x.labels.rotate = FALSE)

# Subset cells / samples
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    cells.use = meta("ident", seurat)!=1)

# Adjust colors
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    colors = c(3,1,2)) #Just changes the color order, probably most useful for dittoDimPlots
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    color.panel = c("red", "orange", "purple"))
```

![](vignettes/QuickStart6_Customizations.png)

```
# Output data
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    data.out = TRUE)

# Add plotly hovering
dittoBarPlot(seurat, "ident", group.by = "RNA_snn_res.0.8",
    do.hover = TRUE)
```

# Color-blindness Friendliness

dittoSeq has many methods to make its plots color-blindness friendly:

### 1. The default color palette is built to work for the most common forms of colorblindness.

I am a protanomalous myself (meaning I am red-green impaired, but more red than green impaired), so I chose colors for dittoSeq that I could tell apart. These colors also work for deuteranomolies (red-green, but more green than red) the most common form of color-blindness.

Note: There are still other forms of colorblindness, tritanomaly (blue deficiency), and complete monochromacy. These are more rare. dittoSeq's default colors are not great for these, but 2 & 3 below can still help!

### 2. Color legend point-sizing is large by default

No color panel can be perfect, but when there are issues, being able to at least establish some of the color differences from the legend helps. For this goal, having the legend examples be large enough is SUPER helpful.

### 3. Lettering overlay

Once the number of colors being used for discrete plotting in `dittoDimPlot` gets too high for even a careful color panel to compensate, letters can be added to by setting `do.letter = TRUE`.

### 4. Shape.by

As an alternate to letting (do.letter & shape.by are incompatible with each other), distinct groups can be displayed using different shapes as well.

### 5. Interactive Plots

Many dittoSeq visualizations offer plotly conversion when a `do.hover` input is set to `TRUE`. Making plots interactive is another great way to make them accessible to individuals with vision impairments. I plan to build such plotly conversion into more functions in the future.

### 6. The **`Simulate`** function

This function allows a cone-typical individual to see what their dittoSeq plot might look like to a colorblind individual.  This function works for all dittoSeq visualizations currently, except for dittoHeatmap.

Note: there are varying degrees of colorblindness. `Simulate` simulates for the most severe cases.

Say this is the code you would use to generate your plot:

```
dittoDimPlot("CD3E", object = seurat, do.letter=F)
```

The code to visualize this as if you were a deuteranope like me is:

```
Simulate(type = "deutan", plot.function=dittoDimPlot, "CD3E", object = seurat, do.letter=F)
```

The Simulate() function's inputs are:

- `type` = "deutan", "protan", "tritan" = the type of colorblindness that you want to simulate.  Deuteranopia is the most common, and involves primarily red color deficiency, and generally also a bit of green.  Protanopia involves primarily green color deficiency, and generally also a bit of red.  Tritanopia involves primarily blue color deficiency.
- `plot.function` = the function you want to use.  R may try to add `()`, but delete that if it does.
- `...` = any and all inputs that go into the plotting function you want to use.

# Demuxlet tools

Included in this package are a set of functions to facilitate Mux-seq applications. For more information on Demuxlet and Mux-sequencing, see the [Demuxlet GitHub Page](https://github.com/statgen/demuxlet). (Impetus: Many Mux-seq experiments will involve generating the side-by-side bulk and single-cell RNAseq data like the rest of the package is built for.)

- **`importDemux()`** - imports Demuxlet info into a pre-made Seurat or SingleCellExperiment object. For more info on its use, see below and `?importDemux` within R.

-	**`demux.calls.summary()`** - Makes a plot of how many calls were made per sample, separated by the separate lanes.  This is very useful for checking the potential accuracy of sample calls when only certain samples went into certain lanes/pools/sequencing runs/etc.  (Note: the default setting is to only show Singlet calls.  Use `singlets.only = FALSE` to include *one of the sample calls* for any doublets.

```
demux.calls.summary(object)
```

-	**`demux.SNP.summary()`** - Useful for checking if you have a lot of cells with very few SNPs. Creates a plot of the number of SNPs per cell that is grouped by individual lane by default.  This function is a simple wrapper for dittoPlot() function with var="demux.N.SNP" and with a number of input defaults adjusted (such as group.by and color.by = "Lane" so that the grouping is done according to 'Lane' metadata.)

```
demux.SNP.summary(object)
```

### `importDemux()` Function:

You will need to point the function to:

- `object` = the target Seurat/SCE object
- `demuxlet.best` = the location(s) of your Demuxlet .best output files.

If your data comes from multiple droplet-gen lanes, then there are two main distinct ways to use the function.

They differ because of specifics of how the data from distinct lanes may have been combined.
See `?importDemux` in R for suggested usage.

#### Metadata created by **`importDemux`**:

Metadata slot name | Description OR the Demuxlet.best column name if directly carried over
--- | ---
Lane | guided by lane.names input, represents of separate droplet-generation lanes, pool, sequencing lane, etc.
Sample | The sample call, from the BEST column
demux.doublet.call | whether the sample was a singlet (SNG), doublet (DBL), or ambiguous (AMB), from the BEST column
demux.RD.TOTL | RD.TOTL
demux.RD.PASS | RD.PASS
demux.RD.UNIQ | RD.UNIQ
demux.N.SNP | N.SNP
demux.PRB.DBL | PRB.DBL
demux.barcode.dup | (Only generated when TRUEs will exist, indicative of a technical issue in the bioinformatics pipeline) whether a cell's barcode referred to only 1 row of the .best file, but multiple distinct cells in the dataset.

#### Summary output:
The import function spits out a quick summary of what was done, which will look something like this:

```
Adding 'Lane' information as meta.data
Extracting the Demuxlet calls
Matching barcodes
Adding Demuxlet info as metadata
Checking for barcode duplicates across lanes...
  No barcode duplicates were found.

SUMMARY:
2 lanes were identified and named:
  Lane1, Lane2
The average number of SNPs per cell for all lanes was: 505.3
Out of 80 cells in the Seurat object, Demuxlet assigned:
    75 cells or 93.8% as singlets
    4 cells or 5% as doublets
    and 1 cells as too ambiguous to call.
0 cells were not annotated in the demuxlet.best file.
```
