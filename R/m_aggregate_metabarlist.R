#' Aggregate PCR replicates
#'
#' Aggregate PCR replicates in a \code{metabarlist} object.
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param groups        a column name in the `pcrs`table corresponding to the sample names to which
#'                      pcr replicates belongs. Default is the `sample_id` column of the table `pcrs` from the
#'                      \code{\link{metabarlist}} object.
#' @param method        a method of aggregation returning a distance matrix. Can be \code{"sum"}, \code{"mean"},
#'                      \code{"prob"}. Default is "sum"
#'
#' @name aggregate_metabarlist
#'
#' @return
#'
#' A \code{metabarlist}
#'
#' @details
#'
#' The function \code{aggregate_metabarlist} is typically used at the end of the data filtration process and aims at aggregating reads and the pcr related information at the sample level. Different methods of aggregation is available:
#'
#' #'\itemize{
#' \item{With method \code{"sum"}, reads of pcr replicates are summed for each OTU}
#' \item{With method \code{"mean"}, reads of pcr replicates are averaged for each OTU. Results are rounded so that to obtain genuine count data}
#' \item{With method \code{"prob"}, the probability of detection is returned for each OTU}
#'
#' For all methods, the `pcr` table is simplified.
#'
#' ## LZ: Fred's coa method to include too?
#' }
#'
#' ##TODO think about other methods of aggregation: e.g. weighted mean, sum and rarefaction, rarefaciton and sum, what else? or let it as a free function (but then this is contradictory with the idea of "simplifying stuff")
#' ##TODO think about how to simplify the pcr table. Wonder if it should not be removed...
#'
#' @examples
#'
#' data(soil_euk)
#'
#' soil_euk_ag <- aggregate_metabarlist(soil_euk, method="sum")
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag)
#'
#' soil_euk_ag <- aggregate_metabarlist(soil_euk, method="prob")
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag) ## output stat are surprising
#'
#' @author Frédéric Boyer, Lucie Zinger
#'
#' @export aggregate_metabarlist


aggregate_metabarlist <- function(metabarlist, groups="sample_id", method="sum") {

  if (suppressWarnings(check_metabarlist(metabarlist))) {

    if (!groups %in% colnames(metabarlist$pcrs)) {
      stop("groups should be in colnames of the metabarlist$pcr table")
    }

    if (!method  %in% c("sum", "mean", "prob")) {
      stop('method should be one of "sum", "mean", or "prob"')
    }

    #aggregate reads
    groups = metabarlist$pcr[,groups]

    if(method == "sum") {
      reads.out <- rowsum(metabarlist$reads, groups)
    } else if(method == "mean") {
      reads.out <- ceiling(rowsum(metabarlist$reads, groups)/as.vector(table(groups)))
    } else if(method == "prob") {
      reads.out <- rowsum(ifelse(metabarlist$reads>0,1,0), groups) / as.vector(table(groups))
    }

    pcr.out <- data.frame(t(sapply(rownames(reads.out), function(x) {
      sub <- metabarlist$pcrs[groups==x,]
      apply(sub, 2, function(x) ifelse(length(unique(x))==1, x[1], NA))
    })))

    pcr.out$plate_no = 1:nrow(pcr.out)
    #not ideal but used to deal with check_metabarlist who requires numerics and non duplicated

    motus.out <- metabarlist$motus[colnames(reads.out),]
    sample.out <- metabarlist$samples[rownames(reads.out[pcr.out$type=="sample",]),]

    out = metabarlist_generator(reads = reads.out, motus = motus.out,
                                pcrs = pcr.out, samples = data.frame(sample.out))

    return(out)
  }
}

