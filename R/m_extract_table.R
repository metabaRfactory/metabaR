#' Extract tables from a \code{metabarlist} object
#'
#' \code{extract_table} is a generic function which extracts a table from a \code{\link{metabarlist}} object.
#'
#'
#' @param metabarlist       a \code{\link{metabarlist}} object
#' @param table             the table required for extraction. Can be one of `reads`, `motus`, `pcrs`, or `samples`
#'
#' @name extract_table
#'
#' @return a numeric matrix `reads` or a dataframe `motus`, `pcrs`, or `samples`
#'
#' @details
#' All \code{metabarlist} objects are composed of four different tables that can be extracted as single tables with this function.
#'
#' @examples
#'
#' data(soil_euk)
#'
#' reads <- extract_table(soil_euk, "reads")
#' motus <- extract_table(soil_euk, "motus")
#'
#' all(colnames(reads) == rownames(motus))
#' is.matrix(reads)
#' is.data.frame(motus)
#' @author Lucie Zinger
#'
#' @export extract_table

extract_table <- function(metabarlist, table = extract_table_methods) {
  if (check_metabarlist(metabarlist)) {
    extract_table_methods <- c("reads", "motus", "pcrs", "samples")
    tab <- match.arg(table, extract_table_methods)

    if (tab == "reads") {
      return(metabarlist$reads)
    } else if (tab == "motus") {
      return(metabarlist$motus)
    } else if (tab == "pcrs") {
      return(metabarlist$pcrs)
    } else {
      return(metabarlist$samples)
    }
  }
}
