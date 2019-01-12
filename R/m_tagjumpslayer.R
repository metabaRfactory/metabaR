#' Filtering potential tag jumps in metabarcoding data
#'
#' Sets to zero abundances of potential tagjumps in a \code{\link{TODEFINE}} object.
#'
#'
#' @param x           a \code{\link{TODEFINE}} object
#' @param threshold   an OTU relative abundance value in an amplicon below which it is considered to be absence
#'
#' @name tagjumpslayer
#'
#' @return a transformed OTU count matrix
#'
#' @details
#' TO WRITE, cite carlsen, esling, etc..
#'
#' @examples
#'
#' data(soil_euk)
#' soil_euk_clean = tagjumpslayer(soil_euk$reads, 0.03)
#' range(colSums(soil_euk$reads))
#' range(colSums(soil_euk_clean))
#' #to improve & add verif plot
#'
#' @author Lucie Zinger
#' @export tagjumpslayer

tagjumpslayer = function(x,threshold=0.03) {

  new = x
  for(y in 1:ncol(x)) {
    cum = cumsum(sort(x[,y,drop=T], decreasing=T))/sum(x[,y])
    if(cum[1]>=(1-threshold)){
      threshold2 = 1-cum[1]
    } else {
      threshold2 = threshold
    }
    out.tmp = 1-cum < threshold2
    out = out.tmp[rownames(x)]
    new[out==T,y] = 0
  }
  return(new)
}
