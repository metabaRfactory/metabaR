#' Filtering potential tag jumps in metabarcoding data
#'
#' Sets to zero abundances of potential tagjumps in a \code{\link{TODEFINE}} object.
#'
#'
#' @param x           a \code{\link{TODEFINE}} object
#' @param threshold   an OTU relative abundance value below which the OTU is considered to be absent
#'
#' @name tagjumpslayer
#'
#' @return a transformed reads count matrix
#'
#' @details
#' Tagjumps are an important bias that lead to the presence of potentially high amounts of false positive in DNA metabarcoding data. The origin of this bias is not well known yet, and currently suspected to be generated during the PCR enrichment process of the sequencing library preparation. Uncomplete PCR amplification at this stage may lead to the formation of chimeras at priming sites, from fragments belonging to two different amplicons. The resulting fragment is therefore strictly identical to the genuine OTU, but its tag combination is artifactual. This bias is also frequence-dependant, i.e. abundant genuine OTUs are more likely to be found in low abundance in samples were they are not supposed to be. The function aims at reducing the amount of such false positives, by considering each OTU separately and set to 0 any abundance representing < 0.03% of the total OTU abundance in the entire dataset.
#'
#' @references Carlsen, T., Aas, A. B., Lindner, D., Vrålstad, T., Schumacher, T., & Kauserud, H. (2012). Don't make a mista (g) ke: is tag switching an overlooked source of error in amplicon pyrosequencing studies?. Fungal Ecology, 5(6), 747-749.
#' @references Esling, P., Lejzerowicz, F., & Pawlowski, J. (2015). Accurate multiplexing and filtering for high-throughput amplicon-sequencing. Nucleic acids research, 43(5), 2513-2524.
#'
#' @references Schnell, I. B., Bohmann, K., & Gilbert, M. T. P. (2015). Tag jumps illuminated–reducing sequence‐to‐sample misidentifications in metabarcoding studies. Molecular ecology resources, 15(6), 1289-1303.
#'
#' @references Zinger, L., Taberlet, P., Schimann, H., Bonin, A., Boyer, F., De Barba, M., ... & Chave, J. (2018). Body size determines soil community assembly in a tropical forest. Molecular ecology.

#' @examples
#'
#' data(soil_euk)
#' soil_euk_clean = tagjumpslayer(soil_euk$reads, 0.03)
#'
#' #identify occurrence of the most abundant OTU
#' idx = which.max(soil_euk$motus$count)
#' p1 = ggpcrplate(attr = soil_euk$reads[,idx],
#'                plate_no = soil_euk$pcrs$plate_no,
#'                plate_col = soil_euk$pcrs$plate_col,
#'                plate_row =  soil_euk$pcrs$plate_row,
#'                control_type = soil_euk$pcrs$Control_type)
#' p1 + scale_size(limits=c(1,max(soil_euk$reads[,idx]))) +
#'      labs(size="# reads", fill="control type") +
#'      scale_fill_manual(values=c("brown", "pink", "cyan4", "red"),
#'                        na.value = "white") +
#'      ggtitle("Distribution of the most abundant OTU")
#'
#'#same on clean data
#'p2 = ggpcrplate(attr = soil_euk_clean[,idx],
#'                plate_no = soil_euk$pcrs$plate_no,
#'                plate_col = soil_euk$pcrs$plate_col,
#'                plate_row =  soil_euk$pcrs$plate_row,
#'                control_type = soil_euk$pcrs$Control_type)
#'p2 + scale_size(limits=c(1,max(soil_euk$reads[,idx]))) +
#'     labs(size="# reads", fill="control type") +
#'     scale_fill_manual(values=c("brown", "pink", "cyan4", "red"),
#'                        na.value = "white") +
#'     ggtitle("Distribution of the most abundant OTU after curation")
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
