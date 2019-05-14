# recursive function to find the non replicating sample
filter_replicat <- function(sub_matrix, threshold) {
  replicat_to_remove <- c()
  if (any(sub_matrix > threshold)) {
    if (nrow(sub_matrix) == 2) {
      replicat_to_remove <- c(replicat_to_remove, rownames(sub_matrix))
    } else {
      replicat_to_remove <- c(colnames(sub_matrix)[which.max(colSums(sub_matrix))],
                              filter_replicat(sub_matrix[-which.max(colSums(sub_matrix)),
                                                         -which.max(colSums(sub_matrix))],
                                              threshold))
    }
  }
  return(replicat_to_remove)
}

# distance function
distance_function <- function(reads) {
  correspondence_analysis <- ade4::dudi.coa(sqrt(reads), scannf=FALSE, nf=2)
  distance_matrix <- dist(correspondence_analysis$li)
  return(distance_matrix)
}

# main function of this script
identify_replicate <- function(metabarlist,
                               FUN = distance_function,
                               groups = metabarlist$pcrs$sample_id,
                               graphics = FALSE) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    if (length(groups) != nrow(metabarlist$pcrs))
      stop('provided groups should have the length of pcrs')

    subset_data <- data.frame(groups=groups, replicating=TRUE,
                             row.names=rownames(metabarlist$pcrs))

    iteration <- 0
    repeat {
      iteration <- iteration+1
      print(paste('Iteration', iteration))

      matrix_with_replicate <- metabarlist$reads[
        rownames(subset_data), ][subset_data$replicating, ]

      distance_matrix <- as.matrix(FUN(matrix_with_replicate))

      replicates <- subset_data[rownames(distance_matrix), 'groups']
      within_replicates <- outer(replicates,
                                 replicates,
                                 FUN = "==") & upper.tri(distance_matrix)
      without_replicates <- outer(replicates,
                                  replicates,
                                  FUN="!=") & upper.tri(distance_matrix)

      if(length(distance_matrix[within_replicates]) < 2){
        stop('Too many replicates have been remove!')
      }
      within_replicate_density <- density(distance_matrix[within_replicates],
                                          from = 0, to = max(distance_matrix),
                                          n = 1000)

      if(length(distance_matrix[without_replicates]) < 2){
        stop('Too many replicates have been remove!')
      }
      without_replicate_density <- density(distance_matrix[without_replicates],
                                           from = 0, to = max(distance_matrix),
                                           n = 1000)

      threshold_distance <- without_replicate_density$x[
        min(which(within_replicate_density$y < without_replicate_density$y))]

      if (graphics) {
        plot(within_replicate_density$x, within_replicate_density$y,
             type = 'l', xlab = 'Distances', ylab = 'Density',
             main = paste('Distances densities\nIteration', iteration))
        lines(without_replicate_density, col = 'blue')
        abline(v = threshold_distance, col = 'red')
      }
      needToBeChecked <- unique(subset_data[
        rownames(which(
          (distance_matrix > threshold_distance) & within_replicates,
          arr.ind = T)),
        'groups'])
      if (length(needToBeChecked) > 0) {
        for (group in needToBeChecked) {
          sub_matrix <-
            distance_matrix[subset_data[rownames(distance_matrix), 'groups'] == group &
                              subset_data[rownames(distance_matrix), 'replicating'],
                            subset_data[rownames(distance_matrix), 'groups'] ==
                              group & subset_data[rownames(distance_matrix), 'replicating']]

          non_replicating <- filter_replicat(sub_matrix,
                                             threshold_distance)
          subset_data$replicating[
            rownames(subset_data) %in% non_replicating] <- FALSE
        }
      }
      else {
        break;
      }
    }
    return(subset_data)
  }
}

# tt <- table(table(subset_data[!subset_data$non_replicating, 'replicates']))
# barplot(tt, main = '#kept replicates')
