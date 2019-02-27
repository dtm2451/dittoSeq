install.packages("devtools")
devtools::install_github("dtm2451/DittoSeq")
library(DittoSeq)
library(Seurat)

#########################################################################################################
# If you downloaded the data linked from my github, this first section was already run. Skip on down.   #
# If not, download from the link below, then run code  below (basically the Satija lab's tutorial code. #
#########################################################################################################

#Data downloaded from:
# https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
# Extracted by double-clicking (Mac is simple for decompressing files =])

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")
#Turn this data into a Seurat object
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,
                           project = "10X_PBMC")
## --- Run "basic" pre-processing
#Create a meta.data slot for the percent reads coming from mitochondrial genes (thought to be
#associated with cell viability at the time of characterization)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
#If you want to view the number of genes per cell, the number of UMI per cell, and the percent.mito, run either of the next lines:
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
multiDBPlot(c("nGene", "nUMI", "percent.mito"), "pbmc", group.by = "orig.ident", color.by = "orig.ident")
#Filter cells based on mitochondiral gene percentage and number of genes captured.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
##Normalization
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000)
##HighlyVariableGeneSelection
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
##Scale (turn into relative values) expression of variable genes
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

## --- Run Dimensional Reduction --- ##

##PCA Calculation
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)
##Run tSNE
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

## --- Run Clustering and identify cell types --- ##

##Find clusters
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10,
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
# -- Note: I'm skipping over a few steps here (namely the markers identification process that went
# how into how the genes these next lines plot were identified as different across clusters
# witihin this dataset.  For full details: https://satijalab.org/seurat/pbmc3k_tutorial.html)
#Some gene exploration from the vignette
#- Seurat data.dimplotting plotting function
FeaturePlot(object = pbmc,
            features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                              "FCGR3A", "LYZ", "PPBP", "CD8A"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")
#- My multiples dimplot plotting function
multiDBDimPlot(c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                 "FCGR3A", "LYZ", "PPBP", "CD8A"),
               "pbmc")
#Cluster renaming:
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells",
                     "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)

#########################
### Data is ready now ###
#########################
# I saved here to generate the file I linked to (or tried to link to if that's broken!)
# Save
# saveRDS(pbmc, file = "~/Downloads/pbmc.rds")

###################################################################################
############ START HERE. If you downloaded the .rds file, START HERE. #############
###################################################################################
# Load #
library(Seurat)
library(DittoSeq)
pbmc <- readRDS("~/Downloads/pbmc.rds")

### What's in this dataset?

#Let's explore the data!!!
# tSNE plot
DBDimPlot("ident", "pbmc", size = 2.5, do.label = T, do.letter = F, label.size = 3)
#Without the labels, with bigger dots, and with letters to help the colorblind
DBDimPlot("ident", "pbmc", size = 2.5, do.letter = NA, label.size = 3) # if do.letter is set to FALSE, feature is turned off.  If do.letter is set to NA, it will turn on when the number of groups is >= 8.

# pca plot
DBDimPlot("ident", "pbmc", reduction.use = "pca") # Because there are 8 groups, lettering actually truns on by default, but we can deactivate that with do.letter =F as I did in one of the lines above.
# With outlines of groups added and lettering off
DBDimPlot("ident", "pbmc", reduction.use = "pca", ellipse = T, do.letter = F)


#Set DEFAULT so we don't have to keep adding "pbmc" / object = "pbmc"
DEFAULT <- "pbmc"
# Now, DBDimPlot requires only 1 variable, DBBarPlot requires 2, and DBPlot requires only 3.


# The idents were assigned celltypes.  Lets make a celltype metadata.
pbmc@meta.data$Celltype <- pbmc@ident

# Let's also make a BroadCelltype metadata, to give us another thing to group things by for playing.
pbmc@meta.data$BroadCelltype <- unlist(sapply(meta("Celltype"), function(X){
  ifelse((X == "B cells" | X =="CD4 T cells" | X =="CD8 T cells" | X =="NK cells"),
         "Lymphoid",
         "Myeloid")
}))

# And finally, let's make 3 mock samples groups
pbmc@meta.data$Sample <- "Sample1"
pbmc@meta.data$Sample[round((1/3*length(pbmc@cell.names)):round((2/3*length(pbmc@cell.names))))] <- "Sample2"
pbmc@meta.data$Sample[round((2/3*length(pbmc@cell.names)+1):length(pbmc@cell.names))] <- "Sample3"
meta.levels("Sample", table.out = T)
#Now the options of metadata are:
get.metas()


# Now we can ask some important questions!
#Do the celltypes express key genes differently?  (These are actually the genes the Satija lab used to call celltypes, so well.. yes)
#- Seurat data.dimplotting plotting function
FeaturePlot(object = pbmc,
            features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                              "FCGR3A", "LYZ", "PPBP", "CD8A"),
            cols.use = c("grey", "blue"),
            reduction.use = "tsne")
######- My multiple dimplot function (in multi form) -#####
multiDBDimPlot(c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                 "FCGR3A", "LYZ", "PPBP", "CD8A"),
               object = "pbmc")
#This line and the one above have the same outcome now!
multiDBDimPlot(c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A",
                 "FCGR3A", "LYZ", "PPBP", "CD8A"))
#We can tweak the range of the colorscale with min/max.  Or change the colors with min.color/max.color.
DBDimPlot("MS4A1", min = 0, max =6, min.color = "grey70", max.color = "purple")

#We can show discrete data, which I did above
DBDimPlot("Celltype")
#We can add labels, highlighted or not, and ellipses, and we can adjust the label size
DBDimPlot("Celltype", do.label = T, do.letter = F, label.size = 2)
DBDimPlot("Celltype", do.label = T, highlight.labels = F, do.letter = F, label.size = 2)
#We can iterate based on values of a metadata, all while still keeping the labels:
multiDBDimPlot_vary_cells("Celltype", cells.use.meta = "BroadCelltype", do.letter = F, do.label = T, label.size = 2)
#The All Cells plot and legend can be turned off, and the number of columns and rows can be adjusted
multiDBDimPlot_vary_cells("Celltype", cells.use.meta = "BroadCelltype", do.letter = F, do.label = T, label.size = 2,
                          all.cells.plot = F, add.single.legend = F,
                          nrow = 1, ncol = 2)  #Note: Default for ncol =3, but Default for nrow = as many as needed
#We can also show continuous data in a plot like these.  It will automatically use the same min/max across all plots.
multiDBDimPlot_vary_cells("MS4A1", cells.use.meta = "BroadCelltype",
                          all.cells.plot = F, add.single.legend = T, legend.title = "MS4A1",
                          nrow = 1, ncol = 3)
#We can add a shape-by variable by setting shape = "metadata"
DBDimPlot("LYZ", shape = "BroadCelltype")
#Or just change all the shapes
DBDimPlot("LYZ", shape = 17) #ggplot uses numbers to refer t different shapes.  For options, see here: https://www.rstudio.com/wp-content/uploads/2015/03/ggplot2-cheatsheet.pdf
#We can also show only certain cells
DBDimPlot("LYZ", cells.use = pbmc@cell.names[meta("BroadCelltype")=="Lymphoid"]) #This can be set with the Seurat method (a list of cell.names)
DBDimPlot("LYZ", cells.use = meta("BroadCelltype")=="Lymphoid") #This can be set with just a logical which would give whether each cell should be kept.  Notice: this is same input that goes into the [] of the cell.names list method

### We can also use plotly and show data about each cell when hovering over them!
library(plotly)
#do.hover = whether to output a ggploty and turn on hover info
#data.hover = a set of gene or metadata names whose data should be shown
DBDimPlot("LYZ", do.hover = T, data.hover = c("LYZ", "Celltype", "BroadCelltype", "Sample"))

#Change the titles and axes labels.
DBDimPlot("LYZ",
          main = "I wanna be the very best",
          sub = "like no one ever was",
          ylab = "to catch them is my real test",
          xlab = "to train them is my cause")
###### All plots: #####
# main = main title   #
# sub = subtitle      #
# ylab = y axis label #
# xlab = x axis label #
#######################

#####- My "y"plotting function (in multi form): (like Seurat's vlnPlot, but A LOT more flexible!) -#####
multiDBPlot(c("MS4A1", "GNLY", "CD3E"),
            object = "pbmc",
            group.by = "Celltype", color.by = "BroadCelltype")
multiDBPlot(c("CD14", "FCER1A", "FCGR3A"),
            object = "pbmc",
            group.by = "Celltype", color.by = "BroadCelltype")
multiDBPlot(c("LYZ", "PPBP", "CD8A"),
            object = "pbmc",
            group.by = "Celltype", color.by = "BroadCelltype")
#We can change how the data is represented:
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("jitter","vlnplot"))  # This is the default
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("vlnplot","jitter"))  # Changing the order sets the plotting order
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("vlnplot","jitter","boxplot"))  # boxplots can also be added
#Boxplots have lots of customizations, which all start with "boxplot."
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("vlnplot","jitter","boxplot"),
       boxplot.width = 0.3,
       boxplot.color = "gray30",
       boxplot.fill = F,
       boxplot.show.outliers = T)  # NOTE: Having outliers on alongside a jitter duplicates those datapoints.
#Jitters have lots of options too, all starting with "jitter.", except for shape.by:
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("vlnplot","jitter"),
       jitter.color = "gray40",
       jitter.size = 1.5,
       shape.by = "Sample",
       jitter.width = 0.3,
       jitter.shape.legend.size = 5)
#We can also reorder, relabel, or rotate the x groupings / labels
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("vlnplot","jitter"),
       rotate.labels = T,
       reorder.x = c(3,5,1,2,7,6,8,4),  # The method: look at an un-reordered plot, then say which position you want the left-most group to move to, then where the originally second-from-left should go,  and so on.
       labels = c("CD4 Ts", "CD8 Ts", "Bs", "NKs", "14 Monos", "FC Monos", "DCs", "Megas"))
#We can also adjust the range and tickmarks of the y-axis
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("jitter","vlnplot"),
       y.breaks = c(0,2.5,5,7.5,10))
#And finally, though there are also some extra tweeks I am leaving out, we can show the raw (or scaled) data as well.
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("jitter","vlnplot"),
       data.type = "raw")
#And we can use do.hover and data.hover to select info to display about the cells if we hover the cursor over the jitter dots
DBPlot("LYZ", group.by = "Celltype", color.by = "BroadCelltype",
       plots = c("jitter","vlnplot"),  # Notice that ggplotly puts the jitter on top regardless of the order of plots.  I may fix this bug in a future version.
       data.type = "raw",
       do.hover = T, data.hover = c("LYZ","CD8A","Celltype","BroadCelltype","Sample"))


#####- My "identity plotting function, DBBarPlot -#######
#### AT THIS TIME, SEURAT DOES NOT HAVE THIS AT ALL #####
DBBarPlot("Celltype", group.by = "Sample") # Look at how there are similar numbers of celltypes in each sample!  (If this weren't because we just randomly split a single sample, it would indicate low inter-individual variation...)
#We can also reorder, relabel, or rotate the x groupings / labels
DBBarPlot("Celltype", group.by = "Sample",
       rotate.labels = F,
       reorder.x = c(3,2,1),  # The method: look at an un-reordered plot, then say which position you want the left-most group to move to, then where the originally second-from-left should go,  and so on.
       x.labels = c("Pikachu", "Eevee", "Miltank"))
#We can adjust what the groupings are called
DBBarPlot("Celltype", group.by = "Sample",
          rename.groups = c("Charizards", "Vaporeons", "Sceptiles", "Zapdos", "Wartortles", "Exeggutors", "Mr. Mimes", "Jynxs"))
#We can also adjust the tickmarks of the y-axis, but NOT the range (always 0 to 1)
DBBarPlot("Celltype", group.by = "Sample",
          y.breaks = seq(0,1,.25))
#We can set a legend title as well as changing all plot titles.
DBBarPlot("Celltype", group.by = "Sample",
          rename.groups = c("Charizards", "Vaporeons", "Sceptiles", "Zapdos", "Wartortles", "Exeggutors", "Mr. Mimes", "Jynxs"),
          legend.title = "Pokemon",
          ylab = "Percent of Wild Encounters",
          main = "Encounters by Visit",
          sub = "Pallet Town Bushes",
          xlab = "Trainers",
          x.labels = c("Ash", "Gary", "Prof Oak"))
#And finally, we can also set do.hover = T to extract the counts/ exact percentage info from the barplot:
DBBarPlot("Celltype", group.by = "Sample",
          rename.groups = c("Charizards", "Vaporeons", "Sceptiles", "Zapdos", "Wartortles", "Exeggutors", "Mr. Mimes", "Jynxs"),
          legend.title = "Pokemon",
          ylab = "Percent of Wild Encounters",
          main = "Encounters by Visit",
          sub = "Pallet Town Bushes",
          xlab = "Trainers",
          x.labels = c("Ash", "Gary", "Prof Oak"),
          do.hover = T) # Notice that ggplotly ignores our rename.groups change of the identity names. I may fix this bug in a future version.

#- My heatmap function, DBHeatmap -#
##### This is a new function that was not part of my original release!!!
library(pheatmap)
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"))

#For real data, you will have more cells than this truncated dataset, so I recommend turning off
# cell clustering when you are trying out tweaks to the look. Do this by adding cluster_cols=FALSE
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cells.annotation = "ident",
          cluster_cols=FALSE)

# Cell names are impossible to read, so I recommend creating a blank meta.data slot, and then setting this to be the cell.names.meta:
pbmc@meta.data$blank <- ""
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank")
#For bulk RNAseq data, or if you have only a few cells, but would like to replace the given names with your samples or a certain metadata slot:
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cells.use = pbmc@cell.names[1:20], #This code directs it to only use the first 20 cells in this dataset
          cell.names.meta = "Sample")

# To zoom in on only a certain set of cells/samples, you can use the cells.use option.  This works the same as for my other functions.
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          cells.annotation = "Celltype",
          cells.use = meta("BroadCelltype")=="Lymphoid")

# Instead of using cell names, `cells.annotation` will annotate cells by cluster or any metadata.  Give this input any "metadata" or "ident"
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          cells.annotation = "ident")
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          cells.annotation = "BroadCelltype")
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          cells.annotation = "Sample")

#To change the colors of the cell annotations, use the input cell.annotation.colors.
#This input requires colors be given in list form, so use list(c("color1","color2","color3"))
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          cells.annotation = "Sample",
          cells.annotation.colors = list(c("orange", "yellow", "purple")))

#To change the colors in the heatmap, you can modify the heatmap.colors input
#Default: heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50)
#This builds a ramp from blue to white to red with 50 steps
#Changing this to heatmap.colors = colorRampPalette(c("black","green"))(30)
#would make the colors go from green to white to darkblue with only 75 steps:
DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
            "FCGR3A","LYZ","PPBP","CD8A"),
          cell.names.meta = "blank",
          heatmap.colors = colorRampPalette(c("green","white","darkblue"))(75))

# If you would like to make any additional tweaks, you can turn on data.out=T and output to a variale in order to output the objects and script.
LIST <- DBHeatmap(c("MS4A1","GNLY","CD3E","CD14","FCER1A",
                    "FCGR3A","LYZ","PPBP","CD8A"),
                  data.out = TRUE)
LIST$pheatmap.script #This is where to find the script


## Other customizations do exist.  Check the documentation for other arguments that are not in here! ##
