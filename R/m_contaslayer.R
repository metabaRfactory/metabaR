#' Detecting contaminants in metabarcoding data using negative controls
#'
#' Uses negative controls to detect contaminant OTUs in a \code{\link{TODEFINE}} object.
#'
#'
#' @param metabarlist a \code{\link{metabarlist}} object
#' @param controls    a vector of amplicon names corresponding to negative controls.
#' @param method      a character string specifying the detection method to be used. Default is \code{"max"}
#'
#' @name contaslayer
#'
#' @return a vector containing the names of OTUs identified as contaminants
#'
#' @details
#' In negative controls, a contaminant should be preferentially amplified since there is no competing DNA. \code{\link{contaslayer}} relies on this assumption and detects OTUs whose relative abundance across the whole dataset is maximum in negative controls.
#' \code{method = "max"} returns the names of OTUs whose frequencies across the entire dataset are maximum in at least one negative control
#' \code{method = "all"} returns the names of OTUs whose frequencies across all negative controls is greater than that across all samples
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # finds contaminants from PCR amplification
#' contaminant <- contaslayer(soil_euk)
#' head(soil_euk$motus[contaminant, ])
#'
#' # Distribution of the most abundant contaminants in the PCR plate design
#' max.conta <- contaminant[which.max(soil_euk$motus[contaminant, "count"])]
#'
#' p <- ggpcrplate(soil_euk,
#'   legend_title = "# reads",
#'   FUN = function(m) {
#'     m$reads[, max.conta]
#'   }
#' )
#' p + scale_size(limits = c(1, max(soil_euk$reads[, max.conta]))) +
#'   ggtitle("Distribution of the most abundant contaminant")
#' @author Lucie Zinger
#' @importFrom vegan decostand
#' @export contaslayer

contaslayer <- function(metabarlist,
                        controls = rownames(metabarlist$pcrs)[which(metabarlist$pcrs$control_type == "pcr")],
                        method = "max") {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads_matrix <- metabarlist$reads
    # transform reads count to frequencies by column (for each reads)
    reads_matrix.fcol <- decostand(reads_matrix, method = "total", MARGIN = 2)

    # get the name of samples having the maximal frequency for each reads id
    reads_matrix.max <- NULL
    for (i in 1:ncol(reads_matrix.fcol)) {
      reads_matrix.max[i] <- rownames(reads_matrix.fcol)[which.max(reads_matrix.fcol[, i])]
    }
    contaminants <- colnames(reads_matrix)[!is.na(match(reads_matrix.max, controls))]

    if (method == "max") {
      return(contaminants)
    } else if (method == "all") {
      idx <- NULL
      for (i in 1:length(contaminants)) {
        reads_id <- contaminants[i]
        controls_sum <- sum(reads_matrix.fcol[controls, reads_id])
        samples_sum <- sum(reads_matrix.fcol[-match(controls, rownames(reads_matrix)), reads_id])
        idx[i] <- controls_sum > samples_sum
      }
      return(contaminants[idx])
    } else {
      stop("The method must be 'max' or 'all'")
    }
  }
}
