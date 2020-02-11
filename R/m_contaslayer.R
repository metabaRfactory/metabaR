#' Separating genuine MOTUs from contaminations using negative controls
#'
#' Uses negative controls to determine if MOTUs in a \code{\link{metabarlist}} object are most likely genuine or most likely contaminants.
#'
#'
#' @param metabarlist   a \code{\link{metabarlist}} object
#' @param method        a character string specifying the detection method to be used. Default: \code{"max"}
#' @param control_types a vector of control types contained in the column control_type of metabarlist$pcrs. This parameter is not used when the parameter \code{controls} is specified. Default is c("pcr", "extraction").
#' @param controls      a vector of amplicon names corresponding to negative controls.
#' @param output_col    a character string for the column name in `metabarlist$motus` in which the result will be stored. Default is "not_contamination"
#'
#' @name contaslayer
#'
#' @return a metabarlist with a new boolean column vector of name `output_col` in `metabarlist$motus` indicating whether MOTUs are genuine (\code{TRUE}) or identified as contaminants (\code{FALSE})
#'
#' @details
#' In negative controls, a contaminant should be preferentially amplified since there is no competing DNA. \code{\link{contaslayer}} relies on this assumption and detects MOTUs whose relative abundance across the whole dataset is maximum in negative controls.
#' \code{method = "max"} considers a MOTU as a contaminant if its frequencies across the entire dataset are maximum in at least one negative control
#' \code{method = "all"} considers a MOTU as a contaminant if its frequencies across all negative controls is greater than that across all samples
#'
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # identify potential genuine from contaminants using PCR or extraction negative controls
#' out <- contaslayer(soil_euk)
#' tail(colnames(out$motus))
#' head(out$motus[which(out$motus$not_contamination==F),])
#' length(which(out$motus$not_contamination==F))
#'
#' # Distribution of the most abundant contaminant MOTU in the PCR plate design
#' max.conta <- rownames(out$motus)[out$motus$not_contamination==F][which.max(out$motus[contaminant, "count"])]
#'
#' p <- ggpcrplate(soil_euk,
#'   legend_title = "# reads",
#'   FUN = function(m) {
#'     m$reads[, max.conta]
#'   }
#' )
#'
#' p + scale_size(limits = c(1, max(soil_euk$reads[, max.conta]))) +
#'   ggtitle("Distribution of the most abundant contaminant")
#'
#'# identify potential genuine from contaminants using extraction negative controls only
#' out <- contaslayer(soil_euk, control_types="extraction", output_col= "not_ext_contamination")
#' tail(colnames(out$motus))
#' length(which(out$motus$not_ext_contamination==F))
#'
#' @author Lucie Zinger, Clement Lionnet
#' @export contaslayer

contaslayer <- function(metabarlist,
                        method = "max",
                        control_types = c("pcr", "extraction"),
                        controls = NULL,
                        output_col = "not_contamination") {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (is.null(controls)) {
      controls <- rownames(metabarlist$pcrs)[which(metabarlist$pcrs$control_type %in% control_types )]
    }
    else if (!all(controls %in% rownames(metabarlist$reads))) {
      v <- unique(controls)
      message(paste(v[!v %in% rownames(metabarlist$reads)],
                    "from 'controls' not found in 'rownames(metabarlist$reads)'",
                    collapse = "\n"
      ))
      stop("All values in 'controls' should have a corresponding entry in 'metabarlist$reads'")
    }

    reads_matrix <- metabarlist$reads
    # transform reads count to frequencies by column (for each reads)
    reads_matrix.fcol <- t(t(reads_matrix)/colSums(reads_matrix))

    # get the name of samples having the maximal frequency for each reads id
    reads_matrix.max <- NULL
    for (i in 1:ncol(reads_matrix.fcol)) {
      reads_matrix.max[i] <- rownames(reads_matrix.fcol)[which.max(reads_matrix.fcol[, i])]
    }
    contaminants <- colnames(reads_matrix)[!is.na(match(reads_matrix.max, controls))]

    if (method == "max") {

      contaminants <- contaminants

    } else if (method == "all") {

      idx <- NULL
      for (i in 1:length(contaminants)) {
        reads_id <- contaminants[i]
        controls_sum <- sum(reads_matrix.fcol[controls, reads_id])
        samples_sum <- sum(reads_matrix.fcol[-match(controls, rownames(reads_matrix)), reads_id])
        idx[i] <- controls_sum > samples_sum
      }

      contaminants <- contaminants[idx]

    } else {
      stop("method must be 'max' or 'all'")
    }

    metabarlist$motus[,output_col] <- !rownames(metabarlist$motus) %in% contaminants

    check_metabarlist(metabarlist)
    return(metabarlist)
  }
}
