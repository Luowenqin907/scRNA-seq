library(devtools)
install_github("immunogenomics/harmony")
library(harmony)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)
library(sctransform)
library(umap)
library(cowplot)
library(harmony)
DefaultAssay(sce1)<-'RNA'
DefaultAssay(sce2)<-'RNA'
sce1$sample<-'sample1'
sce2$sample<-'sample2'
sce<-merge(sce1,sce2,add.cell.ids = c('sample1','sample2'),merge.data=F)
sce.data<-GetAssayData(sce,slot = 'counts')
scenew<- CreateSeuratObject(counts = sce.data)  %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = scenew@var.genes, npcs = 50, verbose = FALSE)
scenew$sample <- sce$sample
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scenew, reduction = "pca", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = scenew, features = "PC_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)

##Run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
scenew <- scenew %>% 
  RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scenew, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = scenew, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = scenew, features = "harmony_1", group.by = "sample", pt.size = .1)
plot_grid(p1,p2)

ElbowPlot(scenew,ndims=50)
scenew <- scenew %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4) %>% 
  identity()

options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(scenew, reduction = "tsne", pt.size = 1.8,label = T,group.by = 'sample')