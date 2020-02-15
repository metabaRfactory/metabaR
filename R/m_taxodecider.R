#' Selecting the good taxonomy
#'
#' Selecting a given taxonomic annotation between those obtained from two different databases
#'
#'
#' @param metabarlist   a \code{\link{metabarlist}} object
#' @param best.db       a vector of two database names, the best one in terms of taxonomic information
#'                      reliability being provided first.
#' @param sim.scores    a vector of two column names in the `motus` table corresponding to the similarity scores
#'                      of each database.
#' @param lineage       a vector of two column names in the `motus` table corresponding to the full taxonomic
#'                      lineage obtained for each database. Should be provided in the same order than `best.db`.
#'                      \code{"centroid"} or \code{"pairwise"}. Default is \code{"centroid"}.
#' @param threshold     a similiarty score threshold above which the annotation of the best database is selected if
#'                      both databases yields high similarity scores.
#'
#' @name taxodecider
#'
#' @details
#'
#' The function \code{taxodecider} allows choosing amongst two taxonomic annotation based on the best identity scores and preferences for a given database (e.g. more reliable taxonomy). All taxonomic information should be stored in `metabarlist$motus`.
#'
#' @return a \code{metabarlist} object with a dataframe motus including the preferred taxonomic assignements
#'
#' @examples
#'
#' \dontrun{
#' data(soil_euk)
#' soil_euk <- silva_annotator(
#'    metabarlist = soil_euk,
#'    silva.path = "~/Documents/workspace/metabaRffe_external_data/lit_euk---ssu---otus.csv",
#'    clust.path =  "~/Documents/workspace/metabaRffe_external_data/lit_euk---ssu---sequence_cluster_map---litiere_euk_cl97_agg_filt.clstr")
#'
#' soil_euk$motus$similarity = soil_euk$motus$similarity/100
#'
#' soil_euk2 <- taxodecider(
#'    metabarlist = soil_euk,
#'    best.db = c("silva", "embl"),
#'    sim.scores = c("similarity", "best_identity.order_filtered_embl_r136_noenv_EUK"),
#'    lineage = c("lineage_silva", "path"),
#'    threshold = 0.9
#' )
#'
#'}
#'
#' @author Anne-Sophie Benoiston, Lucie Zinger
#' @export taxodecider



taxodecider = function(metabarlist, best.db, sim.scores, lineage, threshold){

  if (suppressWarnings(check_metabarlist(metabarlist))) {

  if (!all(sim.scores %in% colnames(metabarlist$motus))) {
    stop("sim.scores should be in colnames of the metabarlist$motus table")
  }

  if (!all(sim.scores %in% colnames(metabarlist$motus))) {
    stop("sim.scores should be in colnames of the metabarlist$motus table")
  }

  x <- metabarlist$motus

  bid <- x[,sim.scores]
  tax <- x[,lineage]
  colnames(bid) <- colnames(tax) <- best.db

  max.db = sapply(1:nrow(bid), function(y){
    dt = data.frame(bid = rank(as.numeric(bid[y,]), na.last=F), # column with the rank of the best identities
                    above = apply(bid[y,], 2, function(z){
                      ifelse(as.numeric(z) > threshold, TRUE, FALSE)}))
    ifelse(sum(dt$above, na.rm = T) == 2, best.db[1], row.names(dt)[which.max(dt$bid)]) # best.db
  })

  idx <- cbind(1:nrow(bid), unname(sapply(max.db, function(y) match(y, colnames(bid)))))
  max.bid <- bid[idx]
  max.tax <- tax[idx]
  metabarlist$motus <- cbind(metabarlist$motus,
                            data.frame(best.db = max.db, best.db.sim = as.numeric(max.bid), best.db.lineage = max.tax))

  check_metabarlist(metabarlist)
  return(metabarlist)
  }
}
