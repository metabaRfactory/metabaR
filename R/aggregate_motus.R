#' Aggregate MOTUs
#'
#' Aggregate MOTUs in a \code{metabarlist} object according to a grouping factor or vector
#'
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param groups        a grouping vector or factor for MOTUs.
#'                      \code{NA} values treated as a group level
#' @param FUN           a function of MOTU aggregation.
#'                      Default is the sum of MOTUs read abundances for each grouping value.
#'
#' @return A \code{metabarlist} where the table `reads` contains MOTUs abundances aggregated according to a grouping vector/factor (e.g. taxonomic assignment at the phylum level), using a method defined in \code{FUN} and where number of columns of tables `reads` will be equal to that of the number of groups in `groups`.
#'
#' @details
#'
#' The function \code{aggregate_motus} is typically used for aggregating MOTUs at a given taxonomic resolution. The user is free to use its own method of aggregation, but the most common aggregation method is to sum reads for each taxa and is therefore pre-encoded in \code{FUN_agg_motus_sum}.
#'
#' After aggregation, the information retained in the `motus` table corresponds to the
#'information of the most abundant MOTU in a given group.
#'
#' @seealso \code{\link{aggregate_pcrs}}, \code{\link{apply}}, \code{\link{aggregate}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' ## With MOTU phylum assignment as the grouping factor and
#' ## default grouping aggregation (sum reads across replicates)
#' soil_euk_ag <- aggregate_motus(soil_euk, groups = soil_euk$motus$phylum_name)
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag)
#'
#' ## With a custom function (here equivalent to FUN_agg_sum,
#' ## i.e. summing all MOTU abundances across groups)
#' soil_euk_ag <- aggregate_motus(soil_euk, groups = soil_euk$motus$phylum_name,
#'                                FUN = function(metabarlist, groups){
#'                                        t(rowsum(t(metabarlist$reads), groups))})
#'
#' @author Lucie Zinger
#'
#' @describeIn aggregate_motus Aggregate MOTUs in a \code{metabarlist} object according to a grouping factor or vector for each pcr.
#' @export aggregate_motus


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
      metabarlist$samples[unique(metabarlist$pcrs$sample_id[pcrs.out$type == "sample"]), , drop=FALSE]

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

#' @describeIn aggregate_motus compute the sum of reads according to a grouping factor or vector for each pcr.
#' @export FUN_agg_motus_sum

FUN_agg_motus_sum <- function(metabarlist, groups) {
  reads.out <- t(rowsum(t(metabarlist$reads), groups))
  return(reads.out)
}



