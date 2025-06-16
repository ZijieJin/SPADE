library(entropy)
library(Seurat)
library(pracma)
library(RcppML)
library(NMF)
library(matrixStats)
library(progress)
library(stringr)
library(parallel)
library(future)
library(future.apply)
library(progressr)
library(shiny)
library(compiler)

options(warn = -1)
options(future.globals.maxSize = 99999999999)
options(show.progress = TRUE)

scCommInit <- function(){
  lr_network <<- read.table(paste0(currentpath, '/scriabin_LR_OmniPath.txt'), sep='\t')
  tfs <<- readRDS(paste0(currentpath, "/dorothea.rds"))
  colnames(lr_network) <- c("from", "to")
  colnames(tfs)[3:4] <- c("from", "to")
  LRpairs_str <<- str_c(lr_network$from, lr_network$to, sep = "---")
  genes <- row.names(expr)
  goodgene <- rep(F, dim(lr_network)[1])
  for (l in 1:length(goodgene)) {
    if (lr_network$from[l] %in% genes & lr_network$to[l] %in% genes) {
      goodgene[l] <- T
    }
  }
  lr_network <<- lr_network[goodgene, ]
  LRpairs_str <<- LRpairs_str[goodgene]
  goodgene <- rep(F, dim(tfs)[1])
  for (l in 1:length(goodgene)) {
    if (tfs$from[l] %in% genes & tfs$to[l] %in% genes) {
      goodgene[l] <- T
    }
  }
  tfs <<- tfs[goodgene & tfs$is_inhibition == 0, ]
}

scComm <- function() {
  scCommInit()
  print("##### Evaluating Weights #####")
  myweights <<- CalcWeights_cmp()
  print("##### Evaluating Cell-Cell Communication Score #####")
  print("##### This may take a long time #####")
  if (numcore > 1){
    zz = list()
    lreach = floor(nrow(lr_network) / numcore)
    for(i in 1:(numcore-1)){
      zz[[i]] = (lreach * (i-1) +1):(lreach * i)
    }
    zz[[numcore]] = (lreach * (numcore-1) + 1):nrow(lr_network)
    
    ccires_part = future_lapply(zz, CalcCCI)
    ccires = ccires_part[[1]]
    for (i in 2:numcore){
      ccires = ccires + ccires_part[[i]]
    }
  }else{
    ccires = CalcCCI_cmp(1:nrow(lr_network))
  }
  return(ccires)
}

CalcWeights <- function() {
  weight1 <- weight3 <- totalweight <- rep(1, dim(lr_network)[1])
  # Calc the entropy of l and r
  print("##### Calculating Weight 1 #####")
  lrgenes = unique(c(lr_network$from, lr_network$to))
  ent = rep(0, length(lrgenes))
  names(ent) = lrgenes
  progress <- progress_bar$new(total = length(ent), clear = F)
  for (i in 1:length(ent)){
    ent[i] = entropy::entropy(as.numeric(expr[lrgenes[i],]))
    progress$tick()
  }
  for (l in 1:dim(lr_network)[1]) {
    weight1[l] <- 1 / (1 + ent[lr_network$from[l]] + ent[lr_network$to[l]])
  }
  weight1[is.na(weight1)] <- 0
  weight1 <- weight1 - min(weight1[weight1 > 0])
  weight1[weight1 < 0] <- 0
  weight1[weight1 > 1] <- 1
  
  # Calc the downstream score
  CalcCor = function(id, expr, tfs){
    tfcor_sub = rep(0, length(id))
    for (l in id) {
      thisv1 <- as.numeric(expr[tfs$from[l], ])
      thisv2 <- as.numeric(expr[tfs$to[l], ])
      tfcor_sub = c(tfcor_sub, cor(thisv1, thisv2))
    }
    return(tfcor_sub)
  }
  print("##### Calculating Weight 2 #####")
  progress <- progress_bar$new(total = dim(tfs)[1], clear = F)
  tfcor <- rep(0, nrow = dim(tfs)[1])
  
  for (l in 1:dim(tfs)[1]) {
    thisv1 <- as.numeric(expr[tfs$from[l], ])
    thisv2 <- as.numeric(expr[tfs$to[l], ])
    tfcor[l] <- cor(thisv1, thisv2)
    progress$tick()
  }
  tfcor[is.na(tfcor)] <- 0
  
  for (l in 1:dim(lr_network)[1]) {
    lcor <- mean(tfcor[which(tfs$to == lr_network$from[l])])
    rcor <- mean(tfcor[which(tfs$to == lr_network$to[l])])
    rcor[is.nan(rcor)] <- 0
    lcor[is.nan(lcor)] <- 0
    weight3[l] <- (lcor + rcor) / 2
  }
  weight3 <- sigmoid(weight3, a = 10)
  totalweight <- weight1 * weight3
  return(list(totalweight = totalweight, weight1 = weight1, weight3 = weight3))
}

CalcCCI <- function(ids) {
  totalweight <- myweights$totalweight
  cciscore <- matrix(0, nrow = dim(expr)[2], ncol = dim(expr)[2])
  rownames(cciscore) <- colnames(cciscore) <- colnames(expr)
  progress <- progress_bar$new(total = length(ids), clear = F)
  exprsd = rowSds(as.matrix(expr))
  names(exprsd) = rownames(expr)
  numcell = ncol(expr)
  for (l in ids) {
    progress$tick()
    thissd <- sqrt(exprsd[lr_network$from[l]] * exprsd[lr_network$to[l]])
    if (thissd < 0.001) {
      next() # 这里thissd不是真正的标准差，卡个界
    }
    print(l)
    exprlij <- sqrt(expr[lr_network$from[l], ] %*% t(expr[lr_network$to[l], ]))
    exprlij_nonzero <- exprlij[exprlij > 0]
    nonzeroprop <- length(exprlij_nonzero) / length(exprlij)
    if (nonzeroprop <= 0.001) {
      next()
    }
    nonzeromean = mean(exprlij_nonzero)
    thismean <- nonzeromean * nonzeroprop
    thissd <- sqrt((mean(exprlij_nonzero^2) - nonzeromean^2 + (1 - nonzeroprop) * nonzeromean^2) * nonzeroprop)
    exprlij <- (exprlij - thismean) / thissd
    exprlij = pmax(exprlij, 0)
    exprlij = pmin(exprlij, 5)
    cciaddmat <- exprlij * totalweight[l]
    cciscore <- cciscore + cciaddmat
  }
  cciscore[cciscore < 0] = 0
  return(cciscore)
}

print('##### Loading Data #####')
args = commandArgs()
currentpath = strsplit(args[4], '=')[[1]][2]
currentpath = str_sub(currentpath, end = max(str_locate_all(currentpath,'/')[[1]][,1]))
if (is.na(currentpath)) currentpath = './'
expr_prefix = str_sub(args[6], end = max(str_locate_all(args[6],'\\.')[[1]][,1]) - 1)

numcore = 1
if (numcore > 1) plan(multisession, workers = numcore)

expr <<- read.csv(args[6], row.names = 1)
genenames = rownames(expr)
genenames = toupper(genenames)
dupli = which(duplicated(genenames))
if (length(dupli) > 0) {
  expr <<- expr[-dupli,]
  genenames = genenames[-dupli]
}
expr <<- as.matrix(expr)
rownames(expr) = genenames

CalcWeights_cmp = cmpfun(CalcWeights)
CalcCCI_cmp = cmpfun(CalcCCI)

res = scComm()
res = round(res * 10000) / 10000
write.csv(res, paste0(expr_prefix, '_CCI.csv'), quote = F)

# currentpath = '/Users/zijie/Desktop/SpatialReconstruction/codes/'
# expr <<- read.csv('/Users/zijie/LXY/data/LIHC/SC_expr.csv', row.names = 1)

