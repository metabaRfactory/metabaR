#' Metabarlist summary
#'
#' summary_metabarlist is a generic function used to produce summary statistics of a \code{\link{metabarlist}} object.
#'
#'
#' @param metabarlist   a \code{\link{metabarlist}} object
#' @param method        type of summary to provided. Should match to "dataset" or "type". If "type", parameters `table` and `index` should be provided. Default is "dataset".
#' @param table         the table where the information on which the aggregation is based. Can be only `motus`, `pcrs`, or `samples`. Default is NULL
#' @param index         a character indicating the name of the element, i.e. columns name on which the aggregation is based. Default is NULL


#'
#' @return Function \code{summary_metabarlist} returns basic summary statistics (nb of elements, reads and MOTUs in total or on average per pcrs) of a \code{\link{metabarlist}} object. The format of the value returned depends on the method used.
#'
#' @details
#' \code{summary_metabarlist} returns basic summary statistics of a \code{\link{metabarlist}} object. The summary returned depends on the `method` used :
#' \itemize{
#' \item{"dataset"}{returns a list of two data.frames: dataset_dimension contains the dimensions of the full metabarlist object, dataset_statistics contains the number of reads, motus in pcrs and samples, as well as average and sd values of reads and motus per samples}
#' \item{"type"}{returns a data.frame similar to the dataset_statistics described above}
#' }
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #dataset summary
#' summary_metabarlist(soil_euk, method = "dataset")
#'
#' #data summary per control type (NA = samples)
#' summary_metabarlist(soil_euk, method = "type", table="pcrs", index="control_type")
#'
#' #data summary per phyla
#' summary_metabarlist(soil_euk, method = "type", table="motus", index="phylum_name")
#'
#' #data summary per Habitat
#' summary_metabarlist(soil_euk, method = "type", table="samples", index="Habitat")
#'
#' @author Lucie Zinger
#' @export summary_metabarlist

summary_metabarlist = function(metabarlist, method="dataset", table=NULL, index=NULL) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    extract_type_methods = c("dataset", "type")
    method = match.arg(method, extract_type_methods)

    if(method == "dataset") {
    dataset_dimension = data.frame(t(sapply(metabarlist, dim)))
    colnames(dataset_dimension) = c("n_row", "n_col")

    idx = which(metabarlist$pcrs$type=="sample")
    size = rowSums(metabarlist$reads)
    rich = rowSums(metabarlist$reads>0)
    dataset_statistics = data.frame(nb_reads = c(sum(metabarlist$reads), sum(metabarlist$reads[idx,])),
                                    nb_motus = c(sum(colSums(metabarlist$reads)>0),
                                                 sum(colSums(metabarlist$reads[idx,])>0)),
                                    avg_reads = c(mean(size),mean(size[idx])),
                                    sd_reads = c(sd(size),sd(size[idx])),
                                    avg_motus = c(mean(rich),mean(rich[idx])),
                                    sd_motus = c(sd(rich),sd(rich[idx])),
                                    row.names = c("pcrs", "samples"))
      return(list(dataset_dimension = dataset_dimension,
                  dataset_statistics = dataset_statistics))
    } else {

      extract_table_methods = c("motus", "pcrs", "samples")
      tab = match.arg(table, extract_table_methods)

      if(length(index)==0 | is.character(index)==F)
        stop("character index of the element to select (i.e. column or row name) should be provided")

      idx = ifelse(is.na(as.vector(metabarlist[[tab]][,index])),
                   "NA", as.vector(metabarlist[[tab]][,index]))

      if(tab == "pcrs" | tab == "motus") {

        size = if(tab=="pcrs") {rowSums(metabarlist$reads)} else {colSums(metabarlist$reads)}
        rich = if(tab=="pcrs") {rowSums(metabarlist$reads>0)} else {colSums(metabarlist$reads>0)}

        dataset_statistics = data.frame(nb_elements = if(tab == "pcrs") {as.vector(table(idx))} else {NA},
                                       nb_reads = as.vector(xtabs(size~idx)),
                                       nb_motus = if(tab == "pcrs") {as.vector(xtabs(rich~idx))
                                         } else {as.vector(table(idx))},
                                       avg_reads = aggregate(size, list(idx), mean)$x,
                                       sd_reads = aggregate(size, list(idx), sd)$x,
                                       avg_motus = aggregate(rich, list(idx), mean)$x,
                                       sd_motus = aggregate(rich, list(idx), sd)$x,
                                       row.names = levels(factor(idx)))
      if(tab == "pcrs") {return(dataset_statistics)} else {return(dataset_statistics[,-1])}

      } else {

        idx1 = match(metabarlist$pcrs$sample_id, rownames(metabarlist$samples), nomatch=0)
        size = rowSums(metabarlist$reads[idx1,])
        rich = rowSums(metabarlist$reads[idx1,]>0)

        dataset_statistics = data.frame(nb_pcrs = as.vector(table(idx[idx1])),
                                       nb_reads = as.vector(xtabs(size~idx[idx1])),
                                       nb_motus = as.vector(xtabs(rich~idx[idx1])),
                                       avg_reads = aggregate(size[idx1], list(idx[idx1]), mean)$x,
                                       sd_reads = aggregate(size[idx1], list(idx[idx1]), sd)$x,
                                       avg_motus = aggregate(rich[idx1], list(idx[idx1]), mean)$x,
                                       sd_motus = aggregate(rich[idx1], list(idx[idx1]), sd)$x,
                                       row.names = levels(factor(idx)))
        return(dataset_statistics)
      }
    }
  }
}




