#' Aggregate PCR replicates
#'
#' Aggregate PCR replicates in a \code{metabarlist} object.
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param replicates    a vector containing the sample names to which each pcr
#'                      replicate belongs and within which they should be aggregated.
#'                      Default is the `sample_id` column from the `pcrs` table.
#' @param FUN           a replicate aggregation function.
#'                      Default is the sum of reads per MOTU across replicates.
#'
#' @name aggregate_pcrs
#'
#' @return A \code{metabarlist} where the table `reads` contains MOTU abundances aggregated according to \code{FUN}. The number of rows of the produced `reads` and `pcrs` tables are equal to that of the `samples` table.
#'
#' @details
#'
#' The function \code{aggregate_pcrs} is typically used at the end of the data filtration process and aims to aggregate reads and the pcr related information at the sample level. The user is free to use their own method of aggregation, but the following are often used and therefore pre-encoded:
#'
#' #'\itemize{
#' \item{\code{"FUN_agg_pcrs_sum"}: the reads of pcr replicates are summed for each MOTU}
#' \item{\code{"FUN_agg_pcrs_mean"}: the reads of pcr replicates are averaged for each MOTU.
#'       Results are rounded so as to obtain genuine count data}
#' \item{\code{"FUN_agg_pcrs_prob"}: the probability of detection is returned for each MOTU.
#'       This method is often used in studies dealing with ancient DNA or diet.}
#' }
#'
#' After aggregation, the information contained in the `pcrs` table is averaged if numeric.
#' If none numeric, information is dereplicated if equal across replicates, or concatenated if not.
#'
#'
#'
#' @seealso \code{\link{aggregate_motus}}, \code{\link{apply}}, \code{\link{aggregate}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' ## With default function (sum reads across replicates)
#' soil_euk_ag <- aggregate_pcrs(soil_euk)
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag)
#'
#' ## With the FUN_agg_prob pre-defined function
#' soil_euk_ag <- aggregate_pcrs(soil_euk, FUN = FUN_agg_pcrs_prob)
#' summary_metabarlist(soil_euk)
#' summary_metabarlist(soil_euk_ag) ## output reads produced do not have much sense in this case.
#'
#' ## With a custom function (here equivalent to FUN_agg_pcrs_sum,
#' ## i.e. summing abundances of all MOTUs across replicates)
#' soil_euk_ag <- aggregate_pcrs(soil_euk,
#'                               FUN = function(metabarlist, replicates){
#'                                  rowsum(metabarlist$reads, replicates)})
#'
#' @author Lucie Zinger, Frédéric Boyer
#'
#' @export aggregate_pcrs
#' @export FUN_agg_pcrs_sum
#' @export FUN_agg_pcrs_mean
#' @export FUN_agg_pcrs_prob


aggregate_pcrs <- function(metabarlist,
                           replicates=NULL,
                           FUN = FUN_agg_pcrs_sum) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {

    if (is.null(replicates)) {
      replicates <- metabarlist$pcrs$sample_id
    }

    if (length(replicates)!=nrow(metabarlist$reads)) {
      stop("`replicates` length should be equal to the number of rows of table `reads`")
    }

    #aggregate reads

    reads.out <- FUN(metabarlist, replicates)

    pcr.out <- data.frame(t(sapply(rownames(reads.out), function(x) {
      sub <- metabarlist$pcrs[replicates == x,]
      apply(sub, 2, function(x) {
        if (is.numeric(x)) {
          mean(x)
        } else {
          ifelse(length(unique(x)) == 1, x[1], paste(unique(sort(x)), collapse = "|"))
        }
      })
    })))

    motus.out <- metabarlist$motus[colnames(reads.out), ]
    sample.out <-
      metabarlist$samples[rownames(reads.out[pcr.out$type == "sample", ]), ]

    metabarlist.out <-
      metabarlist_generator(
        reads = reads.out,
        motus = motus.out,
        pcrs = pcr.out,
        samples = data.frame(sample.out)
      )

    return(metabarlist.out)
  }
}

# sum function
FUN_agg_pcrs_sum <- function(metabarlist, replicates) {
  reads.out <- rowsum(metabarlist$reads, replicates)
  return(reads.out)
}

# mean function
FUN_agg_pcrs_mean <- function(metabarlist, replicates) {
  reads.out <-
    ceiling(rowsum(metabarlist$reads, replicates) / as.vector(table(replicates)))
  return(reads.out)
}

#prob function
FUN_agg_pcrs_prob <- function(metabarlist, replicates) {
  reads.out <-
    rowsum(ifelse(metabarlist$reads > 0, 1, 0), replicates) / as.vector(table(replicates))
  return(reads.out)
}
