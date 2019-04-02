setwd("~/Desktop/Project/data")
library(Seurat)
library(dplyr)
library(hdf5r)
data <- Read10X(data.dir = ".")
dense.size <- object.size(x = as.matrix(x = data))
dense.size#1.08GB
sparse.size <- object.size(x = data)
sparse.size#328.8MB
dense.size/sparse.size#3.4B
pbmc <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

#QC
#nGene: the number of genes for each cell sample
#nUMI: the number of unique transcripts for each cell sample
#percent.mito: the mito transcript percent in each cell sample
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", 'nUMI', "percent.mito"), 
                    low.thresholds = c(2500, 10000, -Inf), high.thresholds = c(8500, 110000, 0.125))

#Normalization
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 4.4, y.cutoff = 0.8)
length(x = pbmc@var.genes)#1380
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

#Dimension Reduction, PCA, Heatmap
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)#Reduce to 20 dimensions
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

#Cluster cells
pbmc2 <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = pbmc2)
pbmc2 <- RunTSNE(object = pbmc2, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = pbmc2)
saveRDS(pbmc2, file = "pbmc.rds")



