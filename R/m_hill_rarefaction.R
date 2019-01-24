#' Generate rarefaction curves using the hill numbers
#'
#' Generate rarefaction curves from a \code{\link{TODEFINE}} object using hill numbers diversity indices, as well as the Good's coverage index.
#'
#'
#' @param x    a community matrix from a \code{\link{TODEFINE}}
#' @param boot nb of resampling to estimate Dq at a given sequencing depth
#' @name hill_rarefaction
#'
#' @return a table of data rarefied at different sequencing depth
#'
#' @details
#' Build a rarefaction analysis of amplicons in a \code{\link{TODEFINE}} object using Hill numbers for q={0,1,2}. These indices are equivalents to richness (q=0), exponential Shannon (q->1), and inverse simpson (q=2). The function also returns the Good's coverage index. But note that this index should be interpreted carefully, as it is based on singletons in each amplicons, some of them being filtered during the bioinformatic process (filtering of absolute singeltons).
#'
#' @references Chao, A., Chiu, C. H., & Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity and differentiation measures through Hill numbers. Annual review of ecology, evolution, and systematics, 45, 297-324.
#'
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #create subset
#' idx = c(grep("-J", soil_euk$pcr$point_id), which(is.na(soil_euk$pcr$Control_type)==F))
#' soil_euk_raref.s = hill_rarefaction(soil_euk$reads[idx,], boot=10)
#'
#' #create a vector differenciating samples of different type matching with \code{hill_rarefaction} output
#' soil_euk$pcr$sample_type = as.vector(soil_euk$pcr$Control_type)
#' soil_euk$pcr$sample_type[is.na(soil_euk$pcr$sample_type)] = "Biological samples"
#' vectype = soil_euk$pcr$sample_type[match(soil_euk_raref.s$sample, rownames(soil_euk$pcr))]
#'
#' gghill_rarefaction(soil_euk_raref.s, vectype)
#'
#' @author Lucie Zinger
#' @importFrom vegan diversity rrarefy
#' @import ggplot2
#' @import reshape2
#' @export hill_rarefaction

hill_rarefaction = function(x, boot = 10) {
  require(vegan)
  out = do.call("rbind", lapply(rownames(x), function(y) {
    #sample sizes
    if(sum(x[y,])>2000) {
      rsize = c(5,10,100,200, 500, 1000, round(seq(2000, sum(x[y,]), length.out = 10)))
    } else {
      rsize = c(1,round(seq(2, sum(x[y,]), length.out = 10)))
    }
    #assess estimates at depths rsize
    out = data.frame(sample = y, reads = rsize, D0 = NA, D1 = NA, D2 = NA, Coverage = NA)
    for(s in 1:length(rsize)) {
      #rarefaction
      r = lapply(1:boot, function(z) {
        rtmp = rrarefy(x[y,], rsize[s])
        rtmp[rtmp>0]})
      out[s,"D0"] = mean(sapply(r, function(z) length(z)))
      out[s, "D1"] = exp(mean(sapply(r, function(z) diversity(z, index = "shannon"))))
      out[s, "D2"] = mean(sapply(r, function(z) diversity(z, index = "invsimpson")))
      out[s, "Coverage"] = mean(sapply(r, function(z) 1-length(which(z==1))/sum(z)))
    }
    return(out)
  }))
  return(out)
}

gghill_rarefaction = function(x, type) {

  b = melt(data.frame(x, type=type), id.var=c("reads", "sample", "type"))
  ggplot(b, aes(x=reads, y=value, group=sample, color=type)) +
    geom_line() +
    facet_wrap(~variable, scale="free", ncol=4) +
    labs(x="#reads", y="estimate", color="Sample type") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "bottom")
}





