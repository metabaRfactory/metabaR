nonReplicatingSample <- function(sub_matrix, threshold) {
  sample_to_remove <- c()
  if (any(sub_matrix > threshold)) {
    if (nrow(sub_matrix) == 2) {
      sample_to_remove <- c(sample_to_remove, rownames(sub_matrix))
    } else {
      sample_to_remove <- c(colnames(sub_matrix)[which.max(colSums(sub_matrix))],
                     nonReplicatingSample(sub_matrix[-which.max(colSums(sub_matrix)),
                      -which.max(colSums(sub_matrix))], threshold))
    }
  }
  return(sample_to_remove)
}


myDistFun <- function(reads, graphics=FALSE, ...) {
  # distFunction
  correspondence_analysis <- ade4::dudi.coa(sqrt(reads), scannf=F, nf=2)
  distance_matrix <- dist(correspondence_analysis$li)
  if (graphics) {
    distM1 <- as.matrix(distance_matrix)
    cols = 1:nlevels(samp[rownames(correspondence_analysis$li), 'replicates'])
    nbReads <- rowSums(dataM[rownames(correspondence_analysis$li), ])
    plot(correspondence_analysis$li, col = cols[as.factor(samp[rownames(correspondence_analysis$li), 'replicates'])],
         pch = 16, cex = nbReads/max(nbReads), ...)
  }
  return(distance_matrix)
}

identify_non_replicating_sample <- function(metabarlist, FUN = myDistFun, graphics = FALSE) {

  samp <- metabarlist$pcrs
  samp$replicates <- samp$sample_id
  samp$nonReplicating <- FALSE

  i <- 0
  repeat {
    i <- i+1
    print(paste('Iteration', i))

    dataM <- metabarlist$reads[rownames(samp), ][! samp$nonReplicating, ]

    distM <- as.matrix(FUN(dataM))


    replicates <- samp[rownames(distM), 'replicates']

    withinReplicates <- outer(replicates, replicates, FUN = "==") & upper.tri(distM)
    notWithinReplicates <- outer(replicates, replicates, FUN="!=") & upper.tri(distM)

    d1 <- density(distM[withinReplicates], from = 0, to = max(distM), n = 1000)
    d2 <- density(distM[notWithinReplicates], from = 0, to = max(distM), n = 1000)

    if (graphics) {
      plot(d1$x, d1$y, type = 'l', xlab = 'Distances', ylab = 'Density',
           main = paste('Distances densities\nIteration', i))
      lines(d2, col = 'red')
      thrDist <- d2$x[min(which(d1$y < d2$y))]
      abline(v = thrDist, col = 'red')
    }

    needToBeChecked <- unique(samp[rownames(which(distM > thrDist & withinReplicates, arr.ind = T)), 'replicates'])
    if (length(needToBeChecked) > 0) {
      for (s in needToBeChecked) {
        pattern <- paste0('^', s)
        sub_matrix <- distM[grep(rownames(distM), pattern = pattern),
                        grep(colnames(distM), pattern = pattern)]
        samp$nonReplicating[rownames(samp) %in% nonReplicatingSample(sub_matrix, thrDist)] <- TRUE
      }
    }
    else {
      break;
    }
  }
  samp[, 'nonReplicating', drop = F]
 }

# tt <- table(table(samp[!samp$nonReplicating, 'replicates']))
# barplot(tt, main = '#kept replicates')
