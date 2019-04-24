#' Subset a \code{metabarlist} object
#'
#' Subsetting a \code{\link{metabarlist}} object
#'
#'
#' @param metabarlist       a \code{\link{metabarlist}} object to be subsetted
#' @param table             the table where the information on which the subsetting is based. Can be only `reads`, `motus`, `pcrs`, or `samples`.
#' @param indices           a vector of names indicating the elements, i.e. rows or columns names to keep in the selected table.
#'
#' @name subset_metabarlist
#'
#' @return a \code{metabarlist} object that contains only the selected elements.
#'
#' @details
#' The subsetting will select particular rows from each table. Factor levels that are unused after selection are dropped.
#' \itemize{
#' \item {If the selection is done on `reads`, `pcrs` or `samples`, the MOTUs not occurring in the selection (i.e. those having a total number of reads of 0) are dropped too.}
#' \item { If the selection is done on `motus`, the pcrs and samples where none of the selected MOTUs are found are kept.}
#'}
#'
#'
#'
#' @examples
#'
#' data(soil_euk)
#'
#'
#' #Create a subset of soil_euk containing only annelids MOTUs
#' ## get motus names assigned to annelids
#' annelids_motus = rownames(soil_euk$motus)[grep("Annelida", soil_euk$motus$path)]
#' ## create the metabarlist object
#' annelids = subset_metabarlist(metabarlist, table = "motus",
#'                               indices = annelids_motus)
#' lapply(annelids, dim)
#'
#' #Create a subset of soil_euk containing only pcrs conducted in plate 1
#' ## get pcrs names from plate 1
#' plate1_pcrs = rownames(soil_euk$pcrs)[which(metabarlist$pcrs$plate_no == 1)]
#' plate1 = subset_metabarlist(metabarlist, table = "pcrs",
#'                             indices = plate1_pcrs)
#'
#' lapply(plate1, dim)
#'
#' #Create a subset of soil_euk containing only positive controls
#' ## get pcrs names corresponding to positive controls
#' poscontrol_pcrs = rownames(soil_euk$pcrs)[which(metabarlist$pcrs$control_type == "positive")]
#' poscontrols = subset_metabarlist(metabarlist, table = "pcrs",
#'                             indices = poscontrol_pcrs)
#'
#' lapply(poscontrols, dim)
#'
#' @author Lucie Zinger
#'
#' @export subset_metabarlist

subset_metabarlist = function(metabarlist, table, indices) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    extract_table_methods = c("reads", "motus", "pcrs", "samples")
    tab = match.arg(table, extract_table_methods)

    if(length(indices)==0 | is.character(indices)==F)
      stop("character indices of the elements to select (i.e. columns or row names) should be provided")

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
        idx = match(as.vector(sort(unique(pcrs$sample_id))), rownames(samples), nomatch=0)
        samples = droplevels(samples[idx,,drop=F])
      } else {
          samples = data.frame(row.names = sort(unique(pcrs$sample_id)))[1:length(unique(pcrs$sample_id)),]
        }

      out = list(reads = reads,
                 motus = motus,
                 pcrs = pcrs,
                 samples = samples)

    } else if(tab == "motus") {

    motus = droplevels(motus[indices,,drop=F])
    reads = reads[,indices,drop=F]

    out = list(reads = reads,
               motus = motus,
               pcrs = pcrs,
               samples = samples)

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


    }

  attr(out, 'class') = "metabarlist"
  check_metabarlist(out)
  if(ncol(samples)==0)
    warning('None of the selected elements are described in the initial samples table;
              an empty sample data.frame is returned')
  return(out)
  }
}
