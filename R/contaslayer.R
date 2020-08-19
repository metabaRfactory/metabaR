#' Detecting contaminants using negative controls
#'
#' Uses negative controls to determine whether MOTUs in a \code{metabarlist} object are more likely to be genuine or contaminants.
#'
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param method        a character string specifying the detection method to be used.
#'                      Default: \code{"max"}
#' @param control_types a vector of control types contained in the column control_type
#'                      of the `pcrs` table. This parameter is not used when the parameter
#'                      \code{controls} is specified.
#'                      Default is c("pcr", "extraction").
#' @param controls      a vector of pcr names corresponding to negative controls.
#' @param output_col    a character string for the column name in the `motus` table, in
#'                      which the result will be stored.
#'                      Default is "not_a_contaminant"
#'
#' @name contaslayer
#'
#' @return a metabarlist with a new boolean column vector of name `output_col` in
#'         the `motus` table indicating whether MOTUs are genuine (\code{TRUE}) or
#'         identified as contaminants (\code{FALSE})
#'
#' @details
#' In negative controls, a contaminant should be preferentially amplified since there is no competing DNA. On the other hand, a MOTU detected in negative controls is not necessarily a contaminant, it can be a genuine MOTU detected in negative controls through tag-jump issues.
#' The function \code{\link{contaslayer}} relies on these assumptions and detects MOTUs whose relative abundance across the whole dataset is highest in negative controls. Two methods are currently available:
#'
#' \itemize{
#' \item{\code{method = "max"} considers a MOTU as a contaminant if its frequencies across the entire dataset are highest in at least one negative control.}
#' \item{\code{method = "all"} considers a MOTU as a contaminant if its frequencies across all negative controls is greater than that across all samples.}
#' }
#'
#' @seealso \code{\link{tagjumpslayer}}, \code{\link{pcrslayer}}
#'          for other data curation procedures.
#'
#' @examples
#'
#' data(soil_euk)
#' library(ggplot2)
#'
#' ## Distinguish genuine MOTUs from contaminants using PCR or extraction negative controls
#' mbl <- contaslayer(soil_euk)
#' tail(colnames(mbl$motus))
#' head(mbl$motus[which(mbl$motus$not_a_contaminant == FALSE),])
#' length(which(mbl$motus$not_a_contaminant == FALSE))
#'
#' ## Distribution of the most abundant contaminant MOTU in the PCR plate design
#' contaminants <- rownames(mbl$motus)[mbl$motus$not_a_contaminant == FALSE]
#' max.conta <- contaminants[which.max(mbl$motus[contaminants, "count"])]
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
#'## Identify contaminants using extraction negative controls only
#' mbl <- contaslayer(soil_euk,
#'                    control_types="extraction",
#'                    output_col= "not_an_ext_contaminant")
#' tail(colnames(mbl$motus))
#' length(which(mbl$motus$not_an_ext_contaminant == FALSE))
#'
#' @author Lucie Zinger, Clement Lionnet
#' @export contaslayer

contaslayer <- function(metabarlist,
                        method = "max",
                        control_types = c("pcr", "extraction"),
                        controls = NULL,
                        output_col = "not_a_contaminant") {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (is.null(controls)) {
      controls <-
        rownames(metabarlist$pcrs)[which(metabarlist$pcrs$control_type %in% control_types)]
    }
    else if (!all(controls %in% rownames(metabarlist$reads))) {
      v <- unique(controls)
      message(
        paste(
          v[!v %in% rownames(metabarlist$reads)],
          "from `controls` not found in the row names of table `reads`",
          collapse = "\n"
        )
      )
      stop("All values in `controls` should be in row names of table `reads`")
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
        samples_sum <-
          sum(reads_matrix.fcol[-match(controls, rownames(reads_matrix)), reads_id])
        idx[i] <- controls_sum > samples_sum
      }

      contaminants <- contaminants[idx]

    } else {
      stop("method must be `max` or `all`")
    }

    metabarlist$motus[,output_col] <- !rownames(metabarlist$motus) %in% contaminants

    check_metabarlist(metabarlist)
    return(metabarlist)
  }
}
