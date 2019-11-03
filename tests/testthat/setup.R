# Tests setup
library(Seurat)
pbmc <- pbmc_small
pbmc.se <- as.SingleCellExperiment(pbmc)
pbmc.rnaseq <- new(
    "RNAseq",
    data = as.matrix(GetAssayData(pbmc)),
    counts = as.matrix(GetAssayData(pbmc, slot = "counts")),
    samples = colnames(pbmc),
    meta.data = pbmc@meta.data)
pbmc.rnaseq <- addDimReduction(
    embeddings = pbmc@reductions$pca@cell.embeddings,
    object = pbmc.rnaseq, name = "pca", key = "RNAseq_PC"
)
