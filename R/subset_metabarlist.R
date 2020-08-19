#' Subsetting metabarlists
#'
#' Return subsets of a \code{\link{metabarlist}} object which meet user-defined conditions.
#'
#'
#' @param metabarlist       a \code{\link{metabarlist}} object
#' @param table             the table where the information on which the subsetting is based.
#'                          Can be one of `reads`, `motus`, `pcrs`, or `samples`.
#' @param indices           a boolean vector indicating the elements,
#'                          i.e. rows or columns to keep.
#'                          Elements to keep should be encoded \code{TRUE}.
#'
#' @name subset_metabarlist
#'
#' @return a \code{\link{metabarlist}} object similar to the initial `metabarlist` except that it contains only the selected elements.
#'
#' @details
#' Subsetting a \code{\link{metabarlist}} will select specific rows from the `reads`, `motus`, `pcrs`, or `samples` tables. Factor levels that are unused after selection are dropped. Note however that:
#' \itemize{
#' \item {If the selection is done on `reads`, `pcrs` or `samples`, the MOTUs not occurring in the pcrs or samples selected (i.e. those having a total number of reads of 0) are dropped.}
#' \item {If the selection is done on `motus`, the pcrs and samples where none of the selected MOTUs are found are kept.}
#' }
#'
#' @seealso \code{\link{droplevels}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#'
#' ## Create a subset of soil_euk containing only annelid MOTUs
#' # i.e. by searching for "Annelida" in MOTUs taxonomic assignments
#' annelids <- subset_metabarlist(soil_euk,
#'                                table = "motus",
#'                                indices = grepl("Annelida", soil_euk$motus$path))
#' summary_metabarlist(annelids)
#'
#' ## Create a subset of soil_euk containing only pcrs conducted in plate 1
#' plate1 <- subset_metabarlist(soil_euk,
#'                              table = "pcrs",
#'                              indices = (soil_euk$pcrs$plate_no == 1))
#' summary_metabarlist(plate1)
#'
#' ## Create a subset of soil_euk containing only positive controls
#' # i.e. excluding biological samples too that are NA in soil_euk$pcrs$control_type
#' poscontrols <- subset_metabarlist(soil_euk,
#'                                   table = "pcrs",
#'                                   indices = (soil_euk$pcrs$control_type == "positive" &
#'                                              !is.na(soil_euk$pcrs$control_type)))
#' summary_metabarlist(poscontrols)
#'
#'
#' @author Lucie Zinger
#'
#' @export subset_metabarlist

subset_metabarlist <- function(metabarlist, table, indices) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    extract_table_methods <- c("reads", "motus", "pcrs", "samples")
    tab <- match.arg(table, extract_table_methods)

    if (length(indices) == 0 |
        is.logical(indices) == F) {
      stop("the elements to select or not (i.e. columns or rows) should be `TRUE/FALSE`")
    }

    reads <- metabarlist$reads
    motus <- metabarlist$motus
    pcrs <- metabarlist$pcrs
    samples <- metabarlist$samples

    if (tab == "reads" | tab == "pcrs") {
      reads <- reads[indices, , drop = F]
      reads <- reads[, colSums(reads) > 0, drop = F]

      motus <- droplevels(motus[colnames(reads), , drop = F])

      pcrs <- droplevels(pcrs[rownames(reads), , drop = F])

      if (any(unique(pcrs$sample_id) %in% rownames(samples))) {
        idx <-
          match(as.vector(sort(unique(pcrs$sample_id))), rownames(samples), nomatch = 0)
        samples <- droplevels(samples[idx, , drop = F])
      } else {
        samples <-
          data.frame(row.names = sort(unique(pcrs$sample_id)))[1:length(unique(pcrs$sample_id)),]
      }

      out <- list(
        reads = reads,
        motus = motus,
        pcrs = pcrs,
        samples = samples
      )
    } else if (tab == "motus") {
      motus <- droplevels(motus[indices, , drop = F])
      reads <- reads[, indices, drop = F]

      out <- list(
        reads = reads,
        motus = motus,
        pcrs = pcrs,
        samples = samples
      )
    } else {
      samples <- droplevels(samples[indices, , drop = F])
      pcrs <- droplevels(pcrs[which(pcrs$sample_id %in% rownames(samples)), , drop = F])
      reads <- reads[rownames(pcrs), , drop = F]
      reads <- reads[, colSums(reads) > 0]
      motus <- droplevels(motus[colnames(reads), , drop = F])

      out <- list(
        reads = reads,
        motus = motus,
        pcrs = pcrs,
        samples = samples
      )
    }

    attr(out, "class") <- "metabarlist"
    check_metabarlist(out)
    if (ncol(samples) == 0) {
      warning(
        "None of the selected elements are described in the initial samples table;
              an empty sample data.frame is returned"
      )
    }
    return(out)
  }
}
