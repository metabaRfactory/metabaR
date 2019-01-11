#' Detects contaminants in metabarcoding data
#'
#' Uses negative controls to detects contaminants OTUs in a \code{\link{TODEFINE}} object.
#'
#'
#' @param x           a \code{\link{TODEFINE}} object
#' @param controls    a vector of amplicon names  corresponding to negative controls.
#' @param method      a character string specifying the detection method to be used. Default is \code{"max"}
#'
#' @name contaslayer
#'
#' @return a vector containing the names of OTUs identified as contaminants
#'
#'#' @details
#' The assumption is that in negative controls, a contaminant should be preferentially amplified as there is no competing DNA.
#' #'  \code{method = "max"} returns OTU names for which frequencies across the entire dataset are maximum in at least one negative control
#'  \code{method = "all"} returns OTU names for which the average frequency across all negative controls is greater than across samples
#'
#'
#'
#' @author Lucie Zinger
#' @export contaslayer

contaslayer = function(x, controls, method="all"){
  require(vegan)
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
      #print(i)
      y = conta[i]
      idx[i] = sum(x.fcol[controls,y]) > sum(x.fcol[-match(controls, rownames(x)),y])
    }
    return(conta[idx])
  }
}
