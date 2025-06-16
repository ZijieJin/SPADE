library(Seurat)
library(ggplot2)
source('/Users/zijie/Desktop/CellCommunication/codes/scComm_v1.2.R')

seobj = readRDS('/Users/zijie/LXY/intestinal/Spatial_A1_adult_with_predictions.RDS')
seobj <- SCTransform(seobj, assay = "Spatial", verbose = FALSE)
seobj <- RunPCA(seobj, assay = "SCT", verbose = FALSE)
seobj <- FindNeighbors(seobj, reduction = "pca", dims = 1:30)
seobj <- FindClusters(seobj, verbose = FALSE)

scdata = readRDS('/Users/zijie/LXY/intestinal/sc.RDS')
scdata$celltype = scdata$Cluster
scdata$cellname = colnames(scdata)
scdata = SetIdent(scdata, value = scdata$celltype)

scdata = subset(scdata, cellname %in% sample(colnames(scdata), 10000))

scdata <- SCTransform(scdata, assay = "RNA", verbose = FALSE)
anchors <- FindTransferAnchors(reference = scdata, query = seobj, normalization.method = 'SCT', recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, 
                                  refdata = scdata$celltype, 
                                  prediction.assay = TRUE,
                                  weight.reduction = seobj[["pca"]], 
                                  dims = 1:30)
zz = predictions.assay$data
majorcelltype = rep('', ncol(seobj))
for (i in 1:length(majorcelltype)){
  majorcelltype[i] = rownames(zz)[which(zz[1:9,i] == zz[10,i])]
}
seobj$spotmajortype = majorcelltype





res = scComm(as.matrix(GetAssayData(seobj, layer = 'counts')), seobj$spotmajortype, 
             lr_database = '/Users/zijie/Desktop/CellCommunication/data/scriabin_LR_OmniPath_extended.txt', 
             tf_database = '/Users/zijie/Desktop/CellCommunication/data/dorothea.rds')

CCImat = res$ccires$cciscore
coord = seobj[["slice1"]]@coordinates
distmat = as.matrix(dist(coord[, c(2,3)]))
correff = rep(0, nrow(distmat))
for (i in 1:nrow(distmat)){
  correff[i] = cor(CCImat[i,], distmat[i,])
}

CCImeasure = rowMeans(CCImat)
summary(correff[CCImeasure > 0.1])

CCImat_norm = scale(CCImat)

df_ori_raw = data.frame(CCI = as.vector(CCImat), DIST = as.vector(distmat))
df_ori = data.frame(CCI = as.vector(CCImat_norm), DIST = as.vector(distmat))

df = df_ori

df_agg = data.frame(DIST = (1:12) * 10, CCI = 0, sd = 0)
for (i in 1:12){
  df_agg$CCI[i] = mean(df_ori_raw$CCI[df_ori_raw$DIST < i * 10 & df_ori_raw$DIST > (i-1) * 10])
  df_agg$sd[i] = sd(df_ori_raw$CCI[df_ori_raw$DIST < i * 10 & df_ori_raw$DIST > (i-1) * 10])
}

df_agg2 = data.frame(CCI = (1:10) / 10, DIST = 0, sd = 0)
for (i in 1:10){
  df_agg2$DIST[i] = mean(df_ori_raw$DIST[df_ori_raw$CCI < i / 10 & df_ori_raw$CCI > (i-1) / 10])
  df_agg2$sd[i] = sd(df_ori_raw$DIST[df_ori_raw$CCI < i / 10 & df_ori_raw$CCI > (i-1) / 10])
}

ggplot(df_agg, aes(x = DIST, y = CCI)) + geom_errorbar(aes(ymin = CCI - 
  sd / 3, ymax = CCI + sd / 3), width = 0.4) + geom_point(size=0.5) + theme_bw() + scale_color_lancet() + 
  theme(axis.text=element_text(size=10, colour = 'black'),axis.title=element_text(size=20,face="bold", colour = 'black'),
        legend.text = element_text(size=10), legend.title = element_text(size=10)) +  labs(x='', y = 'CCI') +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),axis.line = element_line(colour = "black")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(df_agg, aes(x = DIST, y = CCI)) + geom_smooth(se=F) + geom_point(size=0.5) + theme_bw() + scale_color_lancet() + 
  theme(axis.text=element_text(size=10, colour = 'black'),axis.title=element_text(size=20,face="bold", colour = 'black'),
        legend.text = element_text(size=10), legend.title = element_text(size=10)) +  labs(x='', y = 'CCI') +
  theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),axis.line = element_line(colour = "black"))

