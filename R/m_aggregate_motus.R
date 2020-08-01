#' Aggregate MOTUs
#'
#' Aggregate MOTUs in a \code{metabarlist} object according to a grouping factor or vector
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param groups        a grouping vector or factor for MOTUs.
#'                      \code{NA} values treated as a group level
#' @param FUN           a function of MOTU aggregation.
#'                      Default is the sum of MOTUs reads abundances for each grouping value.
#'
#' @name aggregate_motus
#'
#' @return A \code{metabarlist} where the table `reads` contains MOTUs abundances aggregated according to a grouping vector/factor (e.g. taxonomic assignment at the phylum level), using a method defined in \code{FUN} and where number of columns of tables `reads` will be equal to that of the number of groups in `groups`.
#'
#' @details
#'
#' The function \code{aggregate_motus} is typically used for aggregating MOTUs at a given taxonomic resolution. The user is free to use its own method of aggregation, but one method is often used and therefore pre-encoded:
#'
#' #'\itemize{
#' \item{\code{"FUN_agg_motus_sum"}: reads of MOTUs in a given group are summed for each pcr.}
#' }
#'
#' After aggregation, the information contained in the `motus` table corresponds to the
#'information of the most abundant MOTU in a given group.
#'
#' @seealso @seealso \code{\link{aggregate_pcrs}}, \code{\link{apply}}, \code{\link{aggregate}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' ## With MOTU phylum assignment as grouping factor and
#' ## default grouping aggregation (sum reads across replicates)
#' soil_euk_ag <- aggregate_motus(soil_euk)
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag)
#'
#' ## With the FUN_agg_prob pre-defined function
#' soil_euk_ag <- aggregate_motus(soil_euk, FUN = FUN_agg_prob)
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag) ## output on reads do not have much sense in this case.
#'
#' ## With a custom function (here equivalent to FUN_agg_sum,
#' ## i.e. summing all MOTUs abundance across groups)
#' soil_euk_ag <- aggregate_motus(soil_euk,
#'                                FUN = function(reads, groups){
#'                                        t(rowsum(t(metabarlist$reads), groups))})
#'
#' @author Lucie Zinger
#'
#' @export aggregate_motus
#' @export FUN_agg_motus_sum


aggregate_motus <- function(metabarlist,
                            groups=NULL,
                            FUN = FUN_agg_motus_sum) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (is.null(groups)) {
      stop("attribute `groups` is missing and should be provided")
    }

    if (any(is.na(groups))) {
      groups[is.na(groups)] = "NA"
    }

    #aggregate reads

    reads.out <- FUN(metabarlist, groups)

    motus.out <- do.call("rbind", lapply(colnames(reads.out), function(x) {
      sub <- metabarlist$motus[groups == x,]
      sub.count <- colSums(metabarlist$reads[, groups == x, drop=F])
      sub[which.max(sub.count), ]
    }))

    rownames(motus.out) <- colnames(reads.out)

    pcrs.out <- metabarlist$pcrs[rownames(reads.out), ]
    sample.out <-
      metabarlist$samples[rownames(reads.out[pcrs.out$type == "sample", ]), ]

    metabarlist.out <-
      metabarlist_generator(
        reads = reads.out,
        motus = motus.out,
        pcrs = pcrs.out,
        samples = data.frame(sample.out)
      )

    return(metabarlist.out)
  }
}

# sum function
FUN_agg_motus_sum <- function(metabarlist, groups) {
  reads.out <- t(rowsum(t(metabarlist$reads), groups))
  return(reads.out)
}



