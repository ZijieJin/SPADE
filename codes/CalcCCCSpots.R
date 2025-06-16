### 看spot上的CCC
library(Seurat)
library(stringr)

### 一种是直接把spot看成单细胞计算CCC
stCCC = read.csv('/Users/zijie/LXY/CRC_C1/Ypredicted_0_CCI.csv', header = T, row.names = 1)
### 另一种是单细胞根据映射矩阵W把单细胞CCC对应到空间上
stCCC = readRDS('/Users/zijie/LXY/CRC_C1/st_CCI_scW0.rds')

## Read st data
aa = Read10X('/Users/zijie/LXY/CRC_C1/')
image = Read10X_Image('/Users/zijie/LXY/CRC_C1/')
seobj = CreateSeuratObject(aa, assay = 'Spatial')
image <- image[Cells(x = aa)]
DefaultAssay(seobj = image) <- "Spatial"
seobj[["slice1"]] <- image

seobj <- SCTransform(seobj, assay = "Spatial", return.only.var.genes = F)
seobj <- RunPCA(seobj, assay = "SCT", verbose = FALSE)
seobj <- FindNeighbors(seobj, reduction = "pca", dims = 1:30)
seobj <- FindClusters(seobj, verbose = FALSE)
seobj <- RunUMAP(seobj, reduction = "pca", dims = 1:30)

### Calc CCC strength between one spot and nearby spots
distmat = read.csv('/Users/zijie/LXY/CRC_C1/dist.csv', header = T, row.names = 1)
colnames(distmat) = str_replace_all(colnames(distmat), '\\.', '-')
distmat_part = distmat[colnames(seobj), colnames(seobj)]
gaussian_kernel <- exp(- (distmat_part^2) / (2 * 0.03^2))

stCCC[stCCC < 3] = 0
Score = diag(as.matrix(gaussian_kernel) %*% as.matrix(stCCC) + as.matrix(gaussian_kernel) %*% t(as.matrix(stCCC))) / 2
Score = Score / colSums(gaussian_kernel)

seobj = AddMetaData(seobj, data.frame(Score))
SpatialFeaturePlot(seobj, features = 'Score', min.cutoff = 0, max.cutoff = 1)
