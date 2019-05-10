identifyBad <- function(subMat, thr) {
  toSuppress <- c()
  if (any(subMat>thr)) {
    if (nrow(subMat)==2) {
      toSuppress <- c(toSuppress, rownames(subMat))
    } else {
      toSuppress <- c(colnames(subMat)[which.max(colSums(subMat))],
                      identifyBad(subMat[-which.max(colSums(subMat)),
                      -which.max(colSums(subMat))], thr))
    }
  }
  return(toSuppress)
}


samp <- metabarlist$pcrs
samp <- samp[samp$type=="sample",]
samp$replicates <- samp$sample_id
samp$nonReplicating <- FALSE


myDistFun <- function(reads) {
  # distFunction
  h <- ade4::dudi.coa(sqrt(reads), scannf=F, nf=2)
  distM <- as.matrix(dist(h$li))
  #---
  cols = 1:nlevels(samp[rownames(h$li),'replicates'])
  nbReads <- rowSums(dataM[rownames(h$li),])
  plot(h$li, col=cols[as.factor(samp[rownames(h$li),'replicates'])], pch=16,
       main=paste('Correspondance analysis\non', primer,'sqrt transformed data\nIteration',i),
       cex=nbReads/manx(bReads)*10)
  #---
  distM
}

i <- 0
repeat {
  i <- i+1
  print(paste('Iteration',i))

  dataM <- metabarlist$reads[, rownames(samp)][,! samp$nonReplicating]

  distM <- myDistFun(dataM)


  replicates <- samp[rownames(distM), 'replicates']

  withinReplicates <- outer(replicates, replicates, FUN="==") & upper.tri(distM)
  notWithinReplicates <- outer(replicates, replicates, FUN="!=") & upper.tri(distM)

  d1 <- density(distM[withinReplicates], from=0, to=max(distM), n=1000)
  d2 <- density(distM[notWithinReplicates], from=0, to=max(distM), n=1000)

  plot(d1$x, d1$y, type='l', xlab='Distances', ylab='Density',
       main=paste('Distances densities\nIteration',i))
  lines(d2, col='red')
  thrDist <- d2$x[min(which(d1$y<d2$y))]
  abline(v=thrDist, col='red')

  needToBeChecked <- unique(samp[rownames(which(distM>thrDist & withinReplicates, arr.ind=T)), 'replicates'])
  if (length(needToBeChecked)>0) {
    for (s in needToBeChecked) {
      pattern <- paste0('^',s)
      subMat <- distM[grep(rownames(distM), pattern = pattern),
                      grep(colnames(distM), pattern = pattern)]
      samp$nonReplicating[rownames(samp) %in% identifyBad(subMat, thrDist)] <- TRUE
    }
  }
  else {
    break;
  }
}



tt <- table(table(samp[!samp$nonReplicating, 'replicates']))
barplot(tt, main = '#kept replicates')
