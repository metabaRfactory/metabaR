#' Import a ngsfilter file
#'
#' Import a ngsfilter file formatted for the \code{OBITools} package.
#'
#'
#' @param file    path/name of the input
#' @param additional.sep character indicating the field separator in the additional info colum in ngsfilter
#' @param ...    other arguments to be pasted from \code{read.table}
#' @name read_ngsfilter
#'
#' @return a table
#'
#' @details
#' Creates a sample table from a ngsfilter file formatted for the \code{OBITools} package.
#'
#' @examples
#'
#' @author Lucie Zinger

read_ngsfilter = function(file, additional.sep = "=", ...) {
  input = read.table(file, ...)
  colnames(input) = c("experiment", "pcr", "tag_combo", "primerF", "primerR", "additional")

  tags = do.call("rbind", strsplit(as.vector(input$tag_combo), "\\:"))
  colnames(tags) = c("tagF", "tagR")

  additional = gsub("F @ ", "", as.vector(input$additional))
  tmp = strsplit(additional, ";")
  tmp2 = lapply(tmp, function(x) {
    out0 = unlist(strsplit(x, additional.sep))
    out = out0[seq(2, length(out0), by=2)]
    names(out) = out0[seq(1, length(out0), by=2)]
    out
  })
  names.col = unique(sort(unlist(lapply(tmp2, names))))

  d = data.frame()
  for(i in 1:length(tmp2)) for(j in names(tmp2[[i]])) {
    d[i,j] <- tmp2[[i]][j]
  }

  if("position" %in% colnames(d)) {
    d2 = data.frame(plate_no = sapply(strsplit(d$position, "_"), "[[", 1),
                    plate_col = gsub("[A-Z]", "", sapply(strsplit(d$position, "_"), "[[", 2)),
                    plate_row = gsub("[0-9]", "", sapply(strsplit(d$position, "_"), "[[", 2)))
    if(ncol(d)==1) {d = d2} else {d = data.frame(d2, d[,-grep("^position$", colnames(d))])}
  }

  out = data.frame(input[,c("experiment", "pcr", "primerF", "primerR")], tags, d)
  return(out)
}
