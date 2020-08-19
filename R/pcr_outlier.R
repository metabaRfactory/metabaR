#' Identify non-conforming PCR replicates
#'
#' Identify the non-conforming sample or control PCR replicates in the PCRs table from a \code{metabarlist} object.
#' Process numerous reiterations to compare distance densities within PCR replicates and between PCR replicates.
#'
#' @param metabarlist a \code{metabarlist} object
#' @param FUN a function returning a distance matrix. The distance matrix should be a `dist` object which has the same dimensions as the number of rows in the `reads`table, i.e number of PCRs.
#' @param groups a vector containing the replicate identifier. The vector must have the same dimensions as the `PCRs` table from a \code{metabarlist} object. Default = metabarlist$pcrs$sample_id
#' @param graphics a boolean value to plot the distance densities for each iteration. Default = FALSE
#' @param sub_matrix a distance matrix for replicates comparisons
#' @param threshold a threshold below which a pcr is considered as outlier
#'
#' @return a data frame with the replicate groups and a `replicating`column. If not possible, the function will terminate and return an error message.
#'
#' @details
#'
#' This function identifies non-conforming sample or control replicates.
#'
#' The parameter `groups` defines groups of replicates. The vector should be arranged following the format of a `PCRs` table from the \code{metabarlist}.
#' Note: If the distance within replicates is higher than the distance between replicates, the function cannot return any result because all replicates are removed.
#' The parameter `FUN` defines the function used to compute the distance matrix. The function will return an object of class `dist` with the same length as the input table.
#' The default function use the `decostand` and `vegdist` functions from the `vegan` package to perform a correspondance analysis of the `reads` table from the \code{metabarlist}, and returns a distance matrix.
#' Default function detail:
#' bray_function <- function(reads) {
#'   distance_matrix <- vegdist(decostand(reads, method = 'total'), method='bray')
#'   return(distance_matrix)
#' }
#'
#' When the `graphics` parameter is defined as True, a graphic is plotted with the density of distance within replicate and between replicate. The threshold is also plotted as a vertical line at the intersection of the two densities.
#'
#' Note: In the circumstance where numerous sequencing projects have been pooled and analysed on the the same PCR plates, the function must be processed individually for each project, to avoid calculating distances which are meaningless.
#'
#' @examples
#'\dontrun{
#' data(soil_euk)
#'
#' sample_subset <- subset_metabarlist(soil_euk, "pcrs",
#'                                     soil_euk$pcrs$type == "sample")
#' filter_replicat(sample_subset)
#' }
#' @author Frédéric Boyer & Clément Lionnet
#' @describeIn pcr_outlier Identifying the non-replicating samples or controls in the table PCRs from a \code{metabarlist} object.
#' @importFrom stats density
#' @importFrom graphics abline lines
#' @import ade4
#' @import vegan


#' @describeIn pcr_outlier recursive function to find the non replicating samples or controls
filter_replicat <- function(sub_matrix, threshold) {
  replicat_to_remove <- c()
  if (any(sub_matrix > threshold)) {
    if (nrow(sub_matrix) == 2) {
      replicat_to_remove <- c(replicat_to_remove, rownames(sub_matrix))
    } else {
      replicat_to_remove <- c(
        colnames(sub_matrix)[which.max(colSums(sub_matrix))],
        filter_replicat(
          sub_matrix[
            -which.max(colSums(sub_matrix)),
            -which.max(colSums(sub_matrix))
          ],
          threshold
        )
      )
    }
  }
  return(replicat_to_remove)
}

#' @describeIn pcr_outlier distance function with ade4 package and coa analysis
coa_function <- function(reads) {
  correspondence_analysis <- dudi.coa(sqrt(reads), scannf = FALSE, nf = 2)
  distance_matrix <- dist(correspondence_analysis$li)
  return(distance_matrix)
}

#' @describeIn pcr_outlier distance function with vegan package and Bray-Curtis distance
bray_function <- function(reads) {
  distance_matrix <- vegdist(decostand(reads, method = "total"), method = "bray")
  return(distance_matrix)
}

#' @describeIn pcr_outlier main function
pcr_outlier <- function(metabarlist,
                               FUN = bray_function,
                               groups = metabarlist$pcrs$sample_id,
                               graphics = FALSE) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (length(groups) != nrow(metabarlist$pcrs)) {
      stop("provided groups should equal the number of pcrs")
    }

    subset_data <- data.frame(
      groups = groups, replicating = TRUE,
      row.names = rownames(metabarlist$pcrs)
    )

    subset_data[which(rowSums(metabarlist$reads) == 0), "replicating"] <- FALSE

    iteration <- 0
    repeat {
      iteration <- iteration + 1
      print(paste("Iteration", iteration))

      # get only the read for the samples or controls replicating
      matrix_with_replicate <- metabarlist$reads[
        rownames(subset_data),
      ][subset_data$replicating, ]

      # calculate matrix dist
      function_result <- FUN(matrix_with_replicate)

      if (class(function_result) != "dist") {
        stop("The result of provided function is not correct! The function must return object 'dist'!")
      }

      if (length(labels(function_result)) != length(rownames(matrix_with_replicate))) {
        stop("The result of provided function is not correct! The dimension of function result is not correct!")
      }

      if (!all(labels(function_result) %in% rownames(matrix_with_replicate))) {
        stop("The result of provided function is not correct! The labels of function results not correspond to the data!")
      }

      distance_matrix <- as.matrix(function_result)

      replicates <- subset_data[rownames(distance_matrix), "groups"]
      within_replicates <- outer(replicates,
        replicates,
        FUN = "=="
      ) & upper.tri(distance_matrix)
      between_replicates <- outer(replicates,
        replicates,
        FUN = "!="
      ) & upper.tri(distance_matrix)

      if (length(distance_matrix[within_replicates]) < 2) {
        stop("Too many replicates have been remove!")
      }
      within_replicate_density <- density(distance_matrix[within_replicates],
        from = 0, to = max(distance_matrix),
        n = 1000
      )

      if (length(distance_matrix[between_replicates]) < 2) {
        stop("Too many replicates have been remove!")
      }
      between_replicate_density <- density(distance_matrix[between_replicates],
        from = 0, to = max(distance_matrix),
        n = 1000
      )

      #jamais dans les 10 premier %
      threshold_distance <- between_replicate_density$x[min(which(
        cumsum(within_replicate_density$y / sum(within_replicate_density$y)) > 0.1 &
          within_replicate_density$y <= between_replicate_density$y
      ))]


      if (graphics) {
        plot(within_replicate_density$x, within_replicate_density$y,
          type = "l", xlab = "Distances", ylab = "Density",
          main = paste("Distances densities iteration", iteration)
        )
        lines(between_replicate_density, col = "blue")
        abline(v = threshold_distance, col = "red")
      }

      need_to_be_checked <- unique(subset_data[
        rownames(which(
          (distance_matrix > threshold_distance) & within_replicates,
          arr.ind = T
        )),
        "groups"
      ])
      if (length(need_to_be_checked) > 0) {
        for (group in need_to_be_checked) {
          sub_matrix <-
            distance_matrix[
              subset_data[rownames(distance_matrix), "groups"] == group &
                subset_data[rownames(distance_matrix), "replicating"],
              subset_data[rownames(distance_matrix), "groups"] == group &
                subset_data[rownames(distance_matrix), "replicating"]
            ]

          non_replicating <- filter_replicat(
            sub_matrix,
            threshold_distance
          )
          subset_data$replicating[
            rownames(subset_data) %in% non_replicating
          ] <- FALSE
        }
      }
      else {
        break
      }
    }

    #### warning if more than 20% of replicates are removed
    if (dim(subset_data[subset_data$replicating == F, ])[1] / dim(subset_data)[1] > 0.2) {
      warning("More than 20% of replicates are removed !")
    }
    return(subset_data)
  }
}
