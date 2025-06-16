### 选出DE结果中top 100的基因
library(Seurat)
library(stringr)
options(Seurat.object.assay.version = 'v3')

args = commandArgs(trailingOnly = T)


expr = read.csv(args[1], row.names = 1)
allgene = rownames(expr)
expr = as.matrix(expr)
rownames(expr) = allgene
anno = read.csv(args[2])
expr_prefix = str_sub(args[1], end = max(str_locate_all(args[1],'\\.')[[1]][,1]) - 1)

anno$cell_name = str_replace_all(anno$cell_name, '-', '.')

seobj = CreateSeuratObject(expr)
seobj$celltype = anno$cell_type
seobj = SetIdent(seobj, value = seobj$celltype)
m = FindAllMarkers(seobj, logfc.threshold = 0.1)
m1 = m[m$avg_log2FC > 0.15 & m$p_val_adj < 0.05,]
m1$d = 0
for (i in 1:nrow(m1)){
  if (sum(m1$gene == m1$gene[i]) > 5) m1$d[i] = 1
}
m1 = m1[m1$d == 0,]
genelist = data.frame(g = 0, t = 0)
for (t in unique(m1$cluster)){
  mpart = m1[m1$cluster == t,]
  if (nrow(mpart) < 150){
    genelist = rbind(genelist, data.frame(g = mpart$gene, t = t))
  }else{
    q = quantile(mpart$avg_log2FC, 1 - 150 / nrow(mpart))
    genelist = rbind(genelist, data.frame(g = mpart$gene[mpart$avg_log2FC > q], t = t))
  }
}
genelist = genelist[-1,]

write.table(unique(genelist$g), paste0(expr_prefix, '_markers.txt'), sep='\n', quote = F, row.names = F, col.names = F)
