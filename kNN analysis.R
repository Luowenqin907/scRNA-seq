library(Rcpp)
library(RcppArmadillo)
library(scuttle)
library(batchelor)
library(scran)
library(umap)
library(pheatmap)
library(DropletUtils)
library(igraph)
library(SoupX)
library(pbapply)
library(parallel)
library(Seurat)
sce.data<-data.frame(GetAssayData(scenew,slot='counts'))
genes<-rownames(sce.data)
original <- pblapply(unique(scenew$sample), function(t)
{sce_sub <- scenew[, scenew$sample == t]
exprs_in <- as.matrix(GetAssayData(sce_sub))[genes,]
return(exprs_in)})
names(original) <- unique(scenew$sample)
set.seed(100)
mbpca <- multiBatchPCA(scenew$sample1,scenew$sample2)
names(mbpca@listData) <- names(original)
reducedDim(scenew, "PCA") <- do.call(rbind, mbpca)
#similarity assessment using knn regression tool
train = mbpca$sample1
train_celltype <- scenew[, scenew$sample %in% "sample1"]$Celltype
test = mbpca$sample2
test_celltype <-scenew[, scenew$sample %in% "sample2"]$Celltype
k=20

y=as.numeric(train_celltype)
#knn regression...
KNN_regression_output_list <- pblapply(unique(train_celltype), function(c){
  y = train_celltype == c
  y = as.numeric(y)
  KNN_classifier <- FNN::knn.reg(train =as.matrix(train), test =  as.matrix(test), 
                                 k=k, y = y)
  return(KNN_classifier$pred)
})
KNN_regression_output <- data.frame(do.call(cbind, KNN_regression_output_list))
colnames(KNN_regression_output) <- unique(train_celltype)

#aggregate by mean
enrich_data <- aggregate(KNN_regression_output, by =list(factor(test_celltype)), mean)
enrich_data <- data.frame(row.names = enrich_data[, 1], enrich_data[, -1])
pheatmap(enrich_data, color = viridis::magma(20), border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE)