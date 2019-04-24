setwd("~/Desktop/Project/data")
library(Seurat)#v2.3.4
library(dplyr)
library(hdf5r)
data <- Read10X(data.dir = ".")
dense.size <- object.size(x = as.matrix(x = data))
dense.size#1.08GB
sparse.size <- object.size(x = data)
sparse.size#328.8MB
dense.size/sparse.size#3.4B
pbmc <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, project = "BIOINF545")

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
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene","percent.mito"), 
                    low.thresholds = c(-Inf, -Inf), high.thresholds = c(8500, 0.125))

#Normalization
par(mfrow = c(1, 1))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 4.4, y.cutoff = 0.8)##
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          selection.method = dispersion, x.low.cutoff = 0.0125, x.high.cutoff = 4.4, 
                          y.cutoff = 0.5)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          selection.method = "dispersion", top.genes=3000) 
                        
length(x = pbmc@var.genes)#3000
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

#Dimension Reduction, PCA, Heatmap
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
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
pbmc2 <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:20, #1:10
                     resolution = 1, print.output = 0, save.SNN = TRUE)####0.6
PrintFindClustersParams(object = pbmc2)
pbmc2 <- RunTSNE(object = pbmc2, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = pbmc2)
saveRDS(pbmc2, file = "pbmc.rds")

#FeaturePlot
FeaturePlot(object=pbmc, features.plot=c("mCherry"), cols.use=c("pink", "red"), 
            reduction.use="tsne", do.return=T)
FeaturePlot(object=pbmc, features.plot=c("PROM1"), cols.use=c("lightblue", "blue"), 
            reduction.use="tsne", do.return=T)
hist(pbmc@data["mCherry",])
plot(density(pbmc2@data["HIV",]))

# splitting the data and plot
library(ggplot2)
hspc = readRDS('HSPC_pilot_MV.rds')
mid_thresh=3.5
log2_TPM = density(hspc@data["HIV",],from=0,to=max(hspc@data["HIV",]))$x #4294 obs
Density = density(hspc@data["HIV",],from=0,to=max(hspc@data["HIV",]))$y
ViralTranscription = rep(c("Low","High"),length(log2_TPM)/2)
viral_counts = data.frame(log2_TPM,Density,ViralTranscription)
viral_counts$ViralTranscription[log2_TPM <= mid_thresh]="Low"
viral_counts$ViralTranscription[log2_TPM > mid_thresh]="High"

ggplot(viral_counts,aes(x=log2_TPM,y=Density))+geom_line()+
  geom_ribbon(aes(ymin=0,ymax=viral_counts$Density,fill=viral_counts$ViralTranscription))

# really splitting the data
load("/Users/anbyew/Desktop/Project/data/Rdata.RData")
infected <- hspc@data["HIV",] < 3.5 #True: HIV low expression-stem, False: HIV high expression-diff
hspc <- AddMetaData(object = hspc, metadata = infected, col.name = "infected")
hspc.list <- SplitObject(hspc, attribute.1 = "infected")
saveRDS(hspc.list, file = "hspc.list.rds")



#Determine the ‘dimensionality’ of the dataset
hspc <- JackStraw(object = hspc, num.replicate = 100)
hspc <- ScoreJackStraw(object = hspc, dims = 1:20)
JackStrawPlot(object = hspc, dims = 1:15)
ElbowPlot(object = hspc)

hspc <- FindNeighbors(object = hspc, dims = 1:10)
hspc <- FindClusters(object = hspc, resolution = 0.5)
head(x = Idents(object = hspc), 5)
hspc <- RunUMAP(object = hspc, dims = 1:10)
DimPlot(object = hspc, reduction = "umap")

#Cluster cells
hspc <- FindClusters(object = hspc, reduction.type = "pca", dims.use = 1:10, #1:10
                      resolution = 1, print.output = 0, save.SNN = TRUE)####0.6
PrintFindClustersParams(object = hspc)
hspc <- RunTSNE(object = hspc, dims.use = 1:10, do.fast = TRUE)#20
TSNEPlot(object = hspc, do.label=TRUE, label.size=6)

#FeaturePlot
FeaturePlot(object=hspc, features.plot=c("mCherry"), cols.use=c("pink", "red"), 
            reduction.use="tsne", do.return=T)
FeaturePlot(object=hspc, features.plot=c("PROM1"), cols.use=c("lightblue", "blue"), 
            reduction.use="tsne", do.return=T)
FeaturePlot(object=hspc, features.plot=c("HIV"), cols.use=c("yellow", "orange"), 
            +             reduction.use="tsne", do.return=T)


#Add cluster label to metadata, and then split
clusterlb <- c(hspc@ident)
hspc <- AddMetaData(object = hspc, metadata = clusterlb, col.name = "classlb")
hspc.list <- SplitObject(hspc, attribute.1 = "infected")
saveRDS(hspc.list, file = "hspc.list.rds")
saveRDS(hspc, file = "hspc2.rds")


#Split to get the 3000 gene data
hspc3000@data = hspc3000@data[hspc3000@var.genes,]
saveRDS(hspc3000, file = "hspc3000.rds")


#DE!!! #Find DE between (0, 1, 6, 10) VS (rest)
cluster0.markers <- FindMarkers(object = hspc, ident.1 = 0, 
                                ident.2 = c(2,3,4,5,7,8,9,11,12,13,14,15), 
                                min.pct = 0.1)
cluster1.markers <- FindMarkers(object = hspc, ident.1 = 1, 
                                ident.2 = c(2,3,4,5,7,8,9,11,12,13,14,15), 
                                min.pct = 0.1)
cluster6.markers <- FindMarkers(object = hspc, ident.1 = 6, 
                                ident.2 = c(2,3,4,5,7,8,9,11,12,13,14,15), 
                                min.pct = 0.1)
cluster10.markers <- FindMarkers(object = hspc, ident.1 = 10, 
                                ident.2 = c(2,3,4,5,7,8,9,11,12,13,14,15), 
                                min.pct = 0.1)
highlow.markers <- FindMarkers(object = hspc, ident.1 = c(0,1,4,6,10), 
                               ident.2 = c(2,3,5,7,8,9,11,12,13,14,15), 
                                min.pct = 0.25)

hspc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
highlow.markers <- FindMarkers(object = hspc, ident.1 = c(0,1,4,6,10), 
                               ident.2 = c(2,3,5,7,8,9,11,12,13,14,15), 
                               min.pct = 0.25)
highlow.markers[highlow.markers$p_val_adj>0.05,]#180 genes#233

lowhigh.markers <- FindMarkers(object = hspc, ident.1 = c(2,3,5,7,8,9,11,12,13,14,15), 
                               ident.2 = c(0,1,4,6,10), 
                               min.pct = 0.25)
lowhigh.markers[lowhigh.markers$p_val_adj<0.05,]#


VlnPlot(object = hspc, features.plot = c("HOPX", "CD34"))
VlnPlot(object = hspc, features.plot = c("C1QTNF4", "PROM1"))
VlnPlot(object = hspc, features.plot = c("CLDN10", "SPINK2"))

FeaturePlot(object = hspc, features.plot = c("C1QTNF4", "HOPX", "CD34", 
                                             "SPINK2", "PROM1", "CLDN10",
                                             "HIV", "mCherry", "LMNA"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

FeaturePlot(object = hspc, features.plot = c("C1QTNF4", "SPINK2", "EGFL7", "C19orf77",
                                             "HOPX", "CLDN10", "HIV", "mCherry", "LMNA"), 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne")
VlnPlot(object = hspc, features.plot = c("C1QTNF4", "SPINK2", "EGFL7", "C19orf77",
                                         "HOPX", "CLDN10"))

# print out row names for each cluster markers
write.table(row.names(cluster0.markers[cluster0.markers[,"p_val_adj"]<0.05,]), file = "c0.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(row.names(cluster1.markers[cluster1.markers[,"p_val_adj"]<0.05,]), file = "c1.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(row.names(cluster10.markers[cluster10.markers[,"p_val_adj"]<0.05,]), file = "c10.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(row.names(cluster6.markers[cluster6.markers[,"p_val_adj"]<0.05,]), file = "c6.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)

write.table(row.names(highlow.markers[highlow.markers[,"p_val_adj"]<0.05,]), file = "highlow.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)
write.table(row.names(lowhigh.markers[lowhigh.markers[,"p_val_adj"]<0.05,]), file = "lowhigh.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)


highlow.markers['PNMT',]
highlow.markers['HOPX',]
highlow.markers['FAM178B',]
highlow.markers['NEAT1',]
highlow.markers['LGALS1',]
highlow.markers['MYC',]
highlow.markers['MS4A6A',]
highlow.markers['ACY3',]
highlow.markers['MPEG1',]
highlow.markers['CLU',]
highlow.markers['C1QTNF4',]
highlow.markers['SAT1',]
highlow.markers['CKB',]
highlow.markers['HBD',]
highlow.markers['TYROBP',]
highlow.markers['PRSS57',]
highlow.markers['CES1',]
highlow.markers['FABP5',]
highlow.markers['HMGB1',]
highlow.markers['DMTN',]




