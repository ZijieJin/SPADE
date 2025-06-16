library(stringr)
args = commandArgs(trailingOnly = T)
expr_prefix = str_sub(args[1], end = max(str_locate_all(args[1],'\\.')[[1]][,1]) - 1)
coor = read.csv(args[1], row.names = 1)
thisdist = as.matrix(dist(coor))
thisdist = round(thisdist * 100) / 100
rownames(thisdist) = colnames(thisdist) = rownames(coor)
write.csv(thisdist, paste0(expr_prefix, '_dist.csv'), quote = F)