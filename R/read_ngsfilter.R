#' Import an ngsfilter file
#'
#' Import an ngsfilter file formatted for the \code{OBITools} package (Boyer et al. 2016).
#'
#'
#' @param file             path/name of the input ngsfilter txt file
#' @param additional.sep   character indicating the field separator in the additional info column in
#'                         the ngsfilter filter file.
#' @param ...              other arguments to be pasted from \code{read.table}
#' @name read_ngsfilter
#'
#' @return a data.frame
#'
#' @details
#' The function \code{read_ngsfilter} creates a `pcrs` table from an ngsfilter file formatted for the \code{OBITools} package (Boyer et al. 2016).
#'
#' @seealso \code{\link{obifiles_to_metabarlist}}
#'
#' @references Boyer, F., Mercier, C., Bonin, A., Le Bras, Y., Taberlet, P., & Coissac, E. (2016). obitools: a unix-inspired software package for DNA metabarcoding. Molecular Ecology Resources, 16(1), 176-182.
#'
#' @examples
#'
#' \donttest{
#'
#' dir <- tempdir()
#' url = "https://raw.githubusercontent.com/metabaRfactory/metabaR_external_data/master/"
#'
#' ngsfilter_file = "ngsfilter_GWM-768.new_2.txt"
#' ngsfilter_url = paste(url, ngsfilter_file, sep="")
#' ngsfilter_path <- file.path(dir, ngsfilter_file)
#' download.file(ngsfilter_url, ngsfilter_path)
#'
#' ngsfilter <- read_ngsfilter(
#'               file = ngsfilter_path,
#'               sep = "\t",
#'               additional.sep = "=")
#'
#' }
#'
#' @author Lucie Zinger
#' @importFrom utils read.csv2
#' @export read_ngsfilter

read_ngsfilter <- function(file, additional.sep = "=", ...) {
  input <- read.csv2(file, header = F, check.names = F, stringsAsFactors = F, ...)
  colnames(input) <- c("experiment", "pcr_id", "tag_combo", "primer_fwd", "primer_rev", "additional")

  tags <- do.call("rbind", strsplit(as.vector(input$tag_combo), "\\:"))
  colnames(tags) <- c("tag_fwd", "tag_rev")

  additional <- gsub("F @ |T @ ", "", as.vector(input$additional))
  tmp <- strsplit(additional, ";")
  tmp2 <- lapply(tmp, function(x) {
    out0 <- unlist(strsplit(x, additional.sep))
    out <- out0[seq(2, length(out0), by = 2)]
    names(out) <- out0[seq(1, length(out0), by = 2)]
    out
  })
  names.col <- unique(sort(unlist(lapply(tmp2, names))))

  d <- data.frame()
  for (i in 1:length(tmp2)) {
    for (j in names(tmp2[[i]])) {
      d[i, j] <- tmp2[[i]][j]
    }
  }

  if ("position" %in% colnames(d)) {
    d2 <- data.frame(
      plate_no = as.numeric(sapply(strsplit(d$position, "_"), "[[", 1)),
      plate_col = as.numeric(gsub("[A-Z]", "", sapply(strsplit(d$position, "_"), "[[", 2))),
      plate_row = gsub("[0-9]", "", sapply(strsplit(d$position, "_"), "[[", 2))
    )
    if (ncol(d) == 1) {
      d <- d2
    } else {
      d <- data.frame(d2, d[, colnames(d) != "position", drop=F])
    }
  }

  out <- data.frame(input[, c("experiment", "pcr_id", "primer_fwd", "primer_rev")], tags, d)
  return(out)
}
