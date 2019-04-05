#' Subset a \code{metabarlist} object
#'
#' Subsetting a \code{\link{metabarlist}} object
#'
#'
#' @param metabarlist       a \code{\link{metabarlist}} object to be subsetted
#' @param table             the table where the information on which the subsetting is based. Can be only `reads`, `motus`, `pcrs`, or `samples`.
#' @param indices           a numeric vector of indices indicating the elements, rows or columns to keep in the selected table.
#'
#' @name subset_metabarlist
#'
#' @return a \code{metabarlist} object that contains only the selected elements.
#'
#' @details
#' The subsetting will select particular rows from each table. Factor levels that are unused after selection are dropped. If the selection is done on `reads`, `pcrs` or `samples`, the MOTUs not occurring in the selection are dropped too. If the selection is done on `motus`, the pcrs and samples where none of the selected MOTUs are found are kept.
#'
#' @examples
#'
#' data(soil_euk)
#'
#'
#' #create a subset of soil_euk containing only annelids MOTUs
#'
#' annelids = subset_metabarlist(metabarlist, table = "motus",
#'                               indices = grep("Annelida", metabarlist$motus$path))
#' dim(annelids$motus)
#' dim(annelids$reads)
#'
#' #create a subset of soil_euk containing only pcrs conducted in plate 1
#'
#' plate1 = subset_metabarlist(metabarlist, table = "pcrs",
#'                             indices = which(metabarlist$pcrs$plate_no == 1))
#'
#' dim(plate1$reads)
#' dim(plate1$motus)
#'
#' @author Lucie Zinger
#'
#' @export subset_metabarlist

subset_metabarlist = function(metabarlist, table, indices) {

  if(check_metabarlist(metabarlist)) {

    extract_table_methods = c("reads", "motus", "pcrs", "samples")
    tab = match.arg(table, extract_table_methods)

    if(length(indices)==0 | is.numeric(indices)==F)
      stop("numeric indices of the elements to select should be provided")

    reads = metabarlist$reads
    motus = metabarlist$motus
    pcrs = metabarlist$pcrs
    samples = metabarlist$samples

    if(tab == "reads" | tab == "pcrs") {

      reads = reads[indices,,drop=F]
      reads = reads[,colSums(reads)>0,drop=F]

      motus = droplevels(motus[colnames(reads),,drop=F])

      pcrs = droplevels(pcrs[rownames(reads),,drop=F])

      if(any(unique(pcrs$sample_id) %in% rownames(samples))) {
        samples = droplevels(samples[match(as.vector(sort(unique(pcrs$sample_id))),
                                           rownames(samples), nomatch=0),,drop=F])
      } else {
          samples = data.frame(rep(NA, length(unique(pcrs$sample_id))),
                               row.names = sort(unique(pcrs$sample_id)))
        }

      out = list(reads = reads,
                 motus = motus,
                 pcrs = pcrs,
                 samples = samples)
      attr(out, 'class') = "metabarlist"

      if(is.null(samples))
        warning('None of the selected elements are described in the samples table')
      if(any(is.na(samples)))
        warning('Some of the selected elements are not described in the samples table')

      check_metabarlist(out)
      return(out)

    } else if(tab == "motus") {

    motus = droplevels(motus[indices,,drop=F])
    reads = reads[,indices,drop=F]

    out = list(reads = reads,
               motus = motus,
               pcrs = pcrs,
               samples = samples)

    attr(out, 'class') = "metabarlist"
    check_metabarlist(out)
    return(out)

    } else {
      samples = droplevels(samples[indices,,drop=F])
      pcrs = droplevels(pcrs[match(rownames(samples), pcrs$sample_id),,drop=F])
      reads = reads[rownames(pcrs),,drop=F]
      reads = reads[,colSums(reads)>0]
      motus = droplevels(motus[colnames(reads),,drop=F])

      out = list(reads = reads,
                 motus = motus,
                 pcrs = pcrs,
                 samples = samples)

      attr(out, 'class') = "metabarlist"
      check_metabarlist(out)
      return(out)
    }
  }
}
