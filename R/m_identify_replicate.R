#' Identifying the replicate that doesn't replicating
#'
#' Identifying the non-replicating samples or controls in the table PCRs from a \code{\link{metabarlist}} object.
#' Process many iteration to compare the density for distance within replicate and between replicate.
#'
#' @param metabarlist a \code{metabarlist} object
#' @param FUN a function returning a distance matrix. The distance matrix should be a object `dist` which has the same length of table `reads`.
#' @param groups a vector containing the identifier of replicate. The vector must has the same length of the table `PCRs` from a \code{\link{metabarlist}} object. Default = metabarlist$pcrs$sample_id
#' @param graphics a boolean value to plot the distances densities for each iteration. Default = FALSE
#'
#' @name identify_replicate
#'
#' @return data frame with the replicats groups and a column `replicating` or throws a stop
#'
#' @details
#'
#' This function identify the non-replicating samples or controls.
#'
#' The parameter `groups` defined the groups of replicates. The vector should be sorted like the table `PCRs` from the \code{\link{metabarlist}}.
#' Note: if the distance within replicate is too higher than the distance between replicate function cannot return result because all replicates are removed.
#' The parameter `FUN` define the function used to compute the distance matrix. The function must return a object of class `dist` with the same length that its input table.
#' The default function use the function `dudi.coa` of package `ade4` to perform a correspondance analysis of the square root of \code{\link{metabarlist}} table `reads`, and return a matrix distance from the principal components.
#' Default function detail:
#' distance_function <- function(reads) {
#'   correspondence_analysis <- ade4::dudi.coa(sqrt(reads), scannf=FALSE, nf=2)
#'   distance_matrix <- dist(correspondence_analysis$li)
#'   return(distance_matrix)
#' }
#'
#' When the parameter `graphics` is True, a graphic is plotted with the density of distance within replicate and between replicate. The threshold is also plotted like a vertical line at the intersection of two densities.
#'
#' @examples
#'
#' data(soil_euk)
#'
#' sample_subset <- subset_metabarlist(soil_euk, "pcrs", rownames(soil_euk$pcrs)[which(soil_euk$pcrs$type=='sample')])
#' identify_replicate(sample_subset)
#'
#' @author Frédéric Boyer & Clément Lionnet
#' @export identify_replicate


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

      function_result <- FUN(matrix_with_replicate)

      if (class(function_result) != 'dist'){
        stop("The result of provided function is not correct! The function must return object 'dist'!")
      }

      if(length(labels(function_result)) != length(rownames(matrix_with_replicate))){
        stop("The result of provided function is not correct! The dimension of function result is not correct!")
      }

      if(!all(labels(function_result) %in% rownames(matrix_with_replicate))){
        stop("The result of provided function is not correct! The labels of function results not correspond to the data!")
      }

      distance_matrix <- as.matrix(function_result)

      replicates <- subset_data[rownames(distance_matrix), 'groups']
      within_replicates <- outer(replicates,
                                 replicates,
                                 FUN = "==") & upper.tri(distance_matrix)
      between_replicates <- outer(replicates,
                                  replicates,
                                  FUN="!=") & upper.tri(distance_matrix)

      if(length(distance_matrix[within_replicates]) < 2){
        stop('Too many replicates have been remove!')
      }
      within_replicate_density <- density(distance_matrix[within_replicates],
                                          from = 0, to = max(distance_matrix),
                                          n = 1000)

      if(length(distance_matrix[between_replicates]) < 2){
        stop('Too many replicates have been remove!')
      }
      between_replicate_density <- density(distance_matrix[between_replicates],
                                           from = 0, to = max(distance_matrix),
                                           n = 1000)

      threshold_distance <- between_replicate_density$x[
        min(which(within_replicate_density$y < between_replicate_density$y))]

      if (graphics) {
        plot(within_replicate_density$x, within_replicate_density$y,
             type = 'l', xlab = 'Distances', ylab = 'Density',
             main = paste('Distances densities\nIteration', iteration))
        lines(between_replicate_density, col = 'blue')
        abline(v = threshold_distance, col = 'red')
      }

      need_to_be_checked <- unique(subset_data[
        rownames(which(
          (distance_matrix > threshold_distance) & within_replicates,
          arr.ind = T)),
        'groups'])
      if (length(need_to_be_checked) > 0) {
        for (group in need_to_be_checked) {
          sub_matrix <-
            distance_matrix[subset_data[rownames(distance_matrix), 'groups'] == group &
                              subset_data[rownames(distance_matrix), 'replicating'],
                            subset_data[rownames(distance_matrix), 'groups'] == group &
                              subset_data[rownames(distance_matrix), 'replicating']]

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
