library(Seurat)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
stCCC_path <- args[1]
distmat_path <- args[2]
out_path <- args[3]

stCCC = read.csv(stCCC_path, header = TRUE, row.names = 1)
distmat = read.csv(distmat_path, header = TRUE, row.names = 1)
colnames(distmat) = str_replace_all(colnames(distmat), '\\.', '-')
distmat_part = distmat[colnames(seobj), colnames(seobj)]
gaussian_kernel <- exp(- (distmat_part^2) / (2 * 0.03^2))

Score = diag(as.matrix(gaussian_kernel) %*% as.matrix(stCCC) + as.matrix(gaussian_kernel) %*% t(as.matrix(stCCC))) / 2
Score = Score / colSums(gaussian_kernel)

df = data.frame(spotname = rownames(stCCC), Score = Score)
write.csv(df, out_path, quote = F)