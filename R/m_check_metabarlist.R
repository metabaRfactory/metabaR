#' Test that the provided list is a well formed metabaRffe list
#'
#'
#' @param metabarlist a \code{\link{TODEFINE}} object
#'
#' @name check_metabarlist
#'
#' @return TRUE or throws an error
#'
#' @details
#'
#' Check for the properties awaited for a well formed metabaRffe list:
#' \itemize{
#' \item {is a list with three attributes named `reads`, `motus` and `pcrs`}
#' \item {`reads` is a numeric matrix}
#' \item {`motus` and `pcrs` are data.frame}
#' \item {`motus` and `pcrs` have mandatory columns, i.e. tag_fwd, tag_rev, primer_fwd, primer_rev, plate_no, plate_col and plate_row for `pcrs`
#'   and sequence for `motus`}
#' }
#'
#' Moreover, the function issues a warning if any sample or OTU is associated to a count of 0
#'
#' @examples
#'
#' data(soil_euk)
#'
#' check_metabarlist(soil_euk)
#'
#' @author Clément Lionnet & Frédéric Boyer
#' @export check_metabarlist


check_metabarlist <- function(metabarlist) {

  # TODO: add a class attribute (or similar mechanism) to the list to test that it IS a metabaRffe list
  stopifnot(is.list(metabarlist))

  stopifnot(all(c("reads", "motus", "pcrs") %in% names(metabarlist)))

  stopifnot(is.matrix(metabarlist$reads) && is.numeric(metabarlist$reads))

  stopifnot(is.data.frame(metabarlist$motus) && is.data.frame(metabarlist$pcrs))

  stopifnot(all(colnames(metabarlist$reads) %in% rownames(metabarlist$motus)))
  stopifnot(all(rownames(metabarlist$reads) %in% rownames(metabarlist$pcrs)))

  stopifnot('sequence' %in% colnames(metabarlist$motus))

  stopifnot(all(c('tag_fwd', 'tag_rev', 'primer_fwd', 'primer_rev', 'plate_no', 'plate_col', 'plate_row') %in% colnames(metabarlist$pcrs)))

  if (any(rowSums(metabarlist$reads)==0)) {
    warning("Some samples have a count of zero !")
  }
  if (any(colSums(metabarlist$reads)==0)) {
    warning('Some OTUs have a count of zero !')
  }
  TRUE

}
