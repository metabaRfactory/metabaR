#' Detecting contaminants in metabarcoding data using negative controls
#'
#' Uses negative controls to detect contaminants OTUs in a \code{\link{TODEFINE}} object.
#'
#'
#' @param x           a \code{\link{TODEFINE}} object
#' @param controls    a vector of amplicon names corresponding to negative controls.
#' @param method      a character string specifying the detection method to be used. Default is \code{"max"}
#'
#' @name contaslayer
#'
#' @return a vector containing the names of OTUs identified as contaminants
#'
#' @details
#' In negative controls, a contaminant should be preferentially amplified as there is no competing DNA. \code{\link{contaslayer}} relies on this assumption and detects OTUs of which relative abundance across the whole dataset is maximum in negative controls.
#' \code{method = "max"} returns OTU names for which frequencies across the entire dataset are maximum in at least one negative control
#' \code{method = "all"} returns OTU names for which the frequency across all negative controls is greater than that across all samples
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #finds contaminants from PCR amplification
#' pcr.controls = rownames(soil_euk$pcrs)[which(soil_euk$pcrs$Control_type=="PCR")]
#' contaminant = contaslayer(soil_euk$reads, controls = pcr.controls)
#' head(soil_euk$motus[contaminant,])
#'
#' #Distribution of the most abundance contaminant in the PCR plate design
#' max.conta = contaminant[which.max(soil_euk$motus[contaminant, "count"])]
#' p = ggpcrlate(attr = soil_euk$reads[,max.conta], plate_no = soil_euk$pcrs$plate_no,
#'                   plate_col = soil_euk$pcrs$plate_col,
#'                   plate_row =  soil_euk$pcrs$plate_row,
#'                   control_type = soil_euk$pcrs$Control_type)
#' p + scale_size(limits=c(1,max(soil_euk$reads[,max.conta]))) +
#'      labs(size="# reads", fill="control type") +
#'      scale_fill_manual(values=c("brown", "pink", "cyan4", "red"),
#'                        na.value = "white") +
#'      ggtitle("Distribution of the most abundant contaminant")
#'
#' @author Lucie Zinger
#' @importFrom vegan decostand
#' @export contaslayer

contaslayer = function(x, controls, method="max"){

  x.fcol = decostand(x, method='total', MARGIN = 2)

  x.max = NULL
  for (i in 1:ncol(x.fcol)) {
    #print(i)
    x.max[i] = rownames(x.fcol)[which.max(x.fcol[,i])]
  }
  conta = colnames(x)[!is.na(match(x.max,controls))]

  if(method == "max") {
    return(conta)
  } else {
    idx = NULL
    for (i in 1:length(conta)) {
      y = conta[i]
      idx[i] = sum(x.fcol[controls,y]) > sum(x.fcol[-match(controls, rownames(x)),y])
    }
    return(conta[idx])
  }
}

